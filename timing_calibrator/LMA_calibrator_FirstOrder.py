#!/usr/bin/env python3

from random import shuffle

import numpy as np
from scipy.optimize import least_squares
import pyproj

from LoLIM import read_LMA as LMA

from calibrator_tools import delay_fitter

class ITRF_to_local:
    def __init__(self, phase_center, reflatlon):
        self.phase_center = phase_center
        self.reflatlon = reflatlon
        
        RTD = 180.0/np.pi
        
        lat = reflatlon[0]/RTD
        lon = reflatlon[1]/RTD
        self.arg0 = np.array([-np.sin(lon),   -np.sin(lat) * np.cos(lon),   np.cos(lat) * np.cos(lon)])
        self.arg1 = np.array([np.cos(lon) ,   -np.sin(lat) * np.sin(lon),   np.cos(lat) * np.sin(lon)])
        self.arg2 = np.array([         0.0,    np.cos(lat),                 np.sin(lat)])
        
    def __call__(self, itrfpos, out=None):
        
        if out is itrfpos:
            print("out cannot be same as itrfpos in convertITRFToLocal. TODO: make this a real error")
            quit()
            
        if out is None:
            ret = np.empty(itrfpos.shape, dtype=np.double )
        else:
            ret = out
        
        ret[:]  = np.outer(itrfpos[...,0] - self.phase_center[0], self.arg0 )
        ret += np.outer(itrfpos[...,1] - self.phase_center[1], self.arg1 )
        ret += np.outer(itrfpos[...,2] - self.phase_center[2], self.arg2 )
        
        return ret
        
class LatLonAlt_to_LocalXYZ:
    
    def __init__(self, LMA_header, alt_override= None, local=True):
        self.center_LatLonAlt = np.array([ LMA_header.array_lat, LMA_header.array_lon, LMA_header.array_alt ])
        self.local = local
        if alt_override is not None:
            self.center_LatLonAlt[2] = alt_override
            
        self.ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        self.lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        self.transformer = pyproj.Transformer.from_proj( self.lla , self.ecef )
            
        self.center_ITRF = self.geoditic_to_ECEF( self.center_LatLonAlt )
        self.ITRF_CS = ITRF_to_local( self.center_ITRF, self.center_LatLonAlt )
        
    def geoditic_to_ECEF(self, lat_lon_alt):
        # return pyproj.transform(self.lla, self.ecef, lat_lon_alt[1], lat_lon_alt[0], lat_lon_alt[2], radians=False)
        return self.transformer.transform( lat_lon_alt[1], lat_lon_alt[0], lat_lon_alt[2], radians=False)
        
    def __call__(self, LatLonAlt):
        # return LatLonAlt
        
        if len( LatLonAlt.shape ) == 2: ##expect first index is item, second index is lat,lon,alt 
            
            flip = np.swapaxes( LatLonAlt, 0,1)
            out = np.swapaxes( self.geoditic_to_ECEF( flip ), 0,1)
            if self.local:
                out =  self.ITRF_CS(out) #convertITRFToLocal(ITRFs, self.center_ITRF, self.center_LatLonAlt )
            return out
        
        elif len( LatLonAlt.shape ) == 1: ## a single lat lon alt
            out = np.reshape(  self.geoditic_to_ECEF( LatLonAlt ),  (1,3))
            if self.local:
                out = self.ITRF_CS( out ) #convertITRFToLocal(ITRF, self.center_ITRF, self.center_LatLonAlt )
            return out[0]
        
    def LMA_in_bounds(self, source, bounds):
        latlonalt = np.array([ source.latitude, source.longitude, source.altitude ])
        XYZ = self( latlonalt )
        inX = bounds[0][0] <= XYZ[0] <= bounds[0][1]
        inY = bounds[1][0] <= XYZ[1] <= bounds[1][1]
        inZ = bounds[2][0] <= XYZ[2] <= bounds[2][1]
        return inX and inY and inZ 

class LMA_fitter_source:
    def __init__(self, LMA_data, LMA_aux, antennas_to_throw, coord_system):
        
        XYZ = coord_system( np.array([LMA_data.latitude, LMA_data.longitude, LMA_data.altitude]) )
        self.original_XYZT = np.append(XYZ, [ LMA_data.time_of_day ])
        self.original_redChi2 = LMA_data.red_chi_squared

        self.antenna_times = np.array([ t for i,t in enumerate(LMA_aux.peak_times)])
        self.antenna_mask = np.array([ int((T-int(T))*(10**9))>0 for T in self.antenna_times  ], dtype=np.bool)
        
        for i in antennas_to_throw:
            self.antenna_mask[i] = False
            
        self.num_ants = np.sum( self.antenna_mask )
        # self.num_ants = len(self.antenna_times)
     
def errors_to_RMS(errors, params=0):
    E = errors*errors
    S = np.sum( E )
    return np.sqrt( S/(len(errors)-params) )
    
def sort_sources_by_second(sources):
    """this sorts LMA sources into second-aligned bins. Throws out sources that stradel seconds. Also sorts by reduced chi squared"""
    ret = {}
    for s in sources:
        second = int( s.original_XYZT[3] )
        good = True
        for AT, AM in zip(s.antenna_times, s.antenna_mask): ## check if all data is in same second
            if AM:
                ant_sec = int(AT)
                if ant_sec != second:
                    good = False
                    break
                
        if good:
            if second not in ret:
                ret[second] = []
            ret[second].append( s )
   
    all_seconds = sorted( ret.keys() )
    # sources = [ sorted(ret[S], key=lambda s: s.original_redChi2 ) for S in all_seconds ]
    sources = [ ret[S] for S in all_seconds ]
    return np.array(all_seconds), sources
    
    

class calibrator:
    def __init__(self, LMA_fname, max_redChiSquared, min_antennas_per_source, prefered_sources_per_ant, max_sources_per_run, min_sources_per_ant=3, min_antennas_per_run=6, antennas_to_throw=[], verbose=True):
        """
            max_redChiSquared and min_antennas_per_source filter the sources by reduced chi-squared and number of active stations
            antennas_to_throw - names of antennas to throw
            prefered_sources_per_ant - algorithm tries to cluster LMA sources such that each antenna has this many measurments per run
            max_sources_per_run - but will not include more sources then this in a run
            min_sources_per_ant - but, if an antenna has less than min_sources_per_ant in a run, then its delay is assumed to be zero
            min_antennas_per_run - minimum number of antennas to be fitted in run, else that run is skipped
        """
        
        self.LMA_fname = LMA_fname
        self.max_redChiSquared = max_redChiSquared
        self.min_antennas_per_source = min_antennas_per_source
        self.prefered_sources_per_ant = prefered_sources_per_ant
        self.min_sources_per_ant = min_sources_per_ant
        self.antennas_to_throw = antennas_to_throw
        self.min_antennas_per_run = min_antennas_per_run
        self.max_sources_per_run = max_sources_per_run
        self.verbose = verbose
        
        if self.verbose:
            print("reading data")
            
        self.LMA_header, self.LMA_data = LMA.read_LMA_file_data( LMA_file )
        self.AUX_data = self.LMA_header.read_aux_file()
        
        self.coord_system = LatLonAlt_to_LocalXYZ( self.LMA_header, None, local=False )
        
        if self.verbose:
            print("sorting data")
        
        self.antenna_names = [ s.name for s in self.LMA_header.antenna_info_list ]
        self.antenna_indeces_to_throw = [i for i,S in enumerate(self.antenna_names) if S in antennas_to_throw]
    
        antenna_LatLonAlts = np.array( [ [s.lat, s.lon,s.alt] for s in self.LMA_header.antenna_info_list ] )
        self.antenna_localXYZs = self.coord_system( antenna_LatLonAlts )
        
        self.num_antennas = len(self.antenna_localXYZs)
        
        
        LMA_fitters = [ LMA_fitter_source(dat, aux, self.antenna_indeces_to_throw, self.coord_system) for dat,aux in zip(self.LMA_data,self.AUX_data) if dat.red_chi_squared<max_redChiSquared ]
        LMA_fitters = [s for s in LMA_fitters if s.num_ants>=min_antennas_per_source]
            
            
        self.sorted_seconds, self.sources_by_second = sort_sources_by_second( LMA_fitters )
        
        # LMA_fitters.sort(key=lambda s: s.original_redChi2)
        
        self.workspace = delay_fitter( self.antenna_localXYZs, self.max_sources_per_run) ## allocate working memory
        
    def get_number_source_per_second(self):
        return self.sorted_seconds, np.array([ len(S) for S in self.sources_by_second ] )
    
    def process_second(self, second, max_num_runs, ftol=3e-16, gtol=3e-16, max_nfev=1000 ):
        """returns dictionary. Keys are station names. Values are lists of found delays for each run"""
        
        i = np.searchsorted( self.sorted_seconds, second)
        source_list = self.sources_by_second[ i ]
        
        
        ##### select sources for calibration #####
        NumSources_per_run_per_station = np.zeros( (max_num_runs, self.num_antennas), dtype=np.int)
        sources_for_calibration_per_run = [ [] for i in range(num_runs) ]
        
        def could_use_source(source, sources_per_antenna):
            for num,has_ant in zip(sources_per_antenna, source.antenna_mask):
                if has_ant and num < self.prefered_sources_per_ant:
                    return True
            return False
        
        shuffle( source_list )
        # source_list = sorted( source_list, key=lambda s: s.num_ants, reverse=True )
        for s in source_list:
            
            all_good = True
            for SL, NumSources_per_station in zip(sources_for_calibration_per_run,NumSources_per_run_per_station) :
                
                if  len(SL)==self.max_sources_per_run  or  np.all( NumSources_per_station >= self.prefered_sources_per_ant ):
                    ## this run doesn't need more sources
                    continue
                
                all_good = False
                
                if could_use_source(s, NumSources_per_station):
                    NumSources_per_station += s.antenna_mask
                    SL.append( s )
                    break
                
            if all_good:
                break
            
            
            
        #### select reference station ####
        station_goodness = np.zeros( self.num_antennas, dtype=np.int )
        for number_per_station in NumSources_per_run_per_station:
            station_goodness += number_per_station >= self.prefered_sources_per_ant
        referance_antenna_i = np.argmax(station_goodness)
        referance_antenna_name = self.antenna_names[ referance_antenna_i ]
        
        
        
        #### process runs ####
        if self.verbose:
            print('processing second', second, ".", len(source_list), 'sources total')
        
        delay_outputs = {}
        drift_outputs = {}
        for run_i, (SL, numSources_per_station) in enumerate(zip(sources_for_calibration_per_run, NumSources_per_run_per_station)):
            
            ## preamble
            if self.verbose:
                print()
                print()
                print('run', run_i)
            station_can_be_calibrated = numSources_per_station >= self.min_sources_per_ant
            
            if not station_can_be_calibrated[referance_antenna_i]:
                if self.verbose:
                    print('SKIPPING. ref. stat. not included')
                continue
            if np.sum( station_can_be_calibrated ) < self.min_antennas_per_run:
                if self.verbose:
                    print('SKIPPING. not enough stations can be calibrated!')
                continue
            
            
            if self.verbose:
                print(len(SL), "sources")
                print("sources per ant:")
                for sname, num in zip(self.antenna_names, numSources_per_station):
                    print(sname, num)
                print()
            
            
            if len(SL)  > self.max_sources_per_run:
                print( 'ERROR! too many sources! This should never be reached!' )
                quit()
        
            ## set workspace
            antenna_order = [ ant_i for ant_i, (sname, is_good) in enumerate( zip(self.antenna_names,station_can_be_calibrated) ) if sname!=referance_antenna_name and is_good]
            antennaI_to_param_map = np.full(shape=len(self.antenna_names), fill_value=-1, dtype=np.int)
        
            for i,ant_i in enumerate(antenna_order):
                antennaI_to_param_map[ant_i] = i
            
            self.workspace.reset()
            self.workspace.set_delay_indeces( antennaI_to_param_map )
            self.workspace.set_drift_indeces( antennaI_to_param_map )
            self.workspace.set_drift_zero( second )

            for source in SL:
                self.workspace.set_event( source.antenna_times, source.antenna_mask )
            
        
            ## initial fit guess
            number_calibrations = len(antenna_order)
            num_parameters = number_calibrations*2 + 4*len(SL)
            
            initial_guess = np.zeros( num_parameters, dtype=np.double )
            
            i = number_calibrations*2
            for s in SL:
                initial_guess[i:i+4] = s.original_XYZT
                i += 4
                
            if self.verbose:
                initial_residuals = self.workspace.loc_objective_fun( initial_guess[number_calibrations*2:] )
                initial_RMS = errors_to_RMS( initial_residuals, 4*len(SL) )
                print('initial RMS:', initial_RMS)
                
            ## fit location
            ret_loc_only = least_squares( self.workspace.loc_objective_fun, initial_guess[number_calibrations*2:], xtol=None, ftol=ftol, gtol=gtol, x_scale='jac', max_nfev=max_nfev )
            if self.verbose:
                loc_residuals = self.workspace.loc_objective_fun( ret_loc_only.x )
                loc_RMS = errors_to_RMS( loc_residuals, 4*len(SL) )
                print('location-only RMS:', loc_RMS, 'fit msg:', ret_loc_only.message)
            
            
            ## fit delays now
            initial_guess[number_calibrations*2:] = ret_loc_only.x
            ret_delays = least_squares(self.workspace.full_objective_fun, initial_guess[number_calibrations:], xtol=None, ftol=ftol, gtol=gtol, x_scale='jac', max_nfev=max_nfev )
            if self.verbose:
                residuals = self.workspace.full_objective_fun( ret_delays.x )
                RMS = errors_to_RMS( residuals, number_calibrations+4*len(SL) )
                print('delays RMS:', RMS, 'fit msg:', ret_delays.message)
                
                
            ## finally, delays and drifts
            initial_guess[number_calibrations:] = ret_delays.x
            ret = least_squares(self.workspace.firstOrder_objective_fun, initial_guess, xtol=None, ftol=ftol, gtol=gtol, x_scale='jac', max_nfev=max_nfev )
            if self.verbose:
                residuals = self.workspace.firstOrder_objective_fun( ret.x )
                RMS = errors_to_RMS( residuals, number_calibrations*2 + 4*len(SL) )
                print('final RMS:', RMS, 'fit msg:', ret.message)
            
                

            ## place results
            found_drifts = ret.x[:number_calibrations]
            found_delays = ret.x[number_calibrations:number_calibrations*2]
            for ant_i, delay, drift in zip(antenna_order, found_delays, found_drifts):
                antenna = self.antenna_names[ant_i]
                
                if self.verbose:
                    print( antenna, delay, drift )
                    
                if antenna not in delay_outputs:
                    delay_outputs[antenna] = []
                delay_outputs[antenna].append( delay )
                    
                if antenna not in drift_outputs:
                    drift_outputs[antenna] = []
                drift_outputs[antenna].append( drift )
                
        return delay_outputs, drift_outputs
        

if __name__ == "__main__":
    # LMA_file = "/home/brian/processed_files/2018/D20180813T153001.413Z/LMA/LOFAR_180813_152000_0600.dat.gz"
    # min_antennas_per_source = 7
    # max_redChiSquared = 0.25
    
    
    LMA_file = "./DATA/COLMA_180522_030000_0600.dat.gz"
    # LMA_file = "./COLMA_180518_221000_0600_cal3.dat.gz"
    min_antennas_per_source = 9
    max_redChiSquared = 1*( (25/70)**2 )  ## NOTE: if this is too small, then only sources whose random error cancels out station error will be used. DO NOT SET TOO SMALL!
    
    sources_per_ant = 30
    max_sources_per_run = 30
    min_sources_per_ant = 1
    min_sources_in_second = 100
    
    num_runs = 10
    max_nfev = 1000
    
    out_fname = './results/COLMA_180522_030000_0600_DriftsCal.txt'
    
    # mode = 'test' ## use this mode to test callibrator on one second of data
    mode = 'full_run' ## use this mode to run callibrator for multiple seconds and ouput result to a file
    
    CAL = calibrator(LMA_fname = LMA_file, 
                     max_redChiSquared = max_redChiSquared, 
                     min_antennas_per_source = min_antennas_per_source, 
                     prefered_sources_per_ant  = sources_per_ant, 
                     max_sources_per_run = max_sources_per_run, 
                     min_sources_per_ant = min_sources_per_ant, 
                     min_antennas_per_run=6, 
                     antennas_to_throw=[], 
                     verbose=True)
    
    all_seconds, number_sources = CAL.get_number_source_per_second()
    
    
    ### RUN ###
    
    if mode == 'full_run':
        out_file = open(out_fname, 'w')
        
        out_file.write(str( len(all_seconds) ))
        out_file.write( ' ' )
        for ant_name in CAL.antenna_names:
            out_file.write( ant_name )
            out_file.write( ' ' )
        out_file.write( '\n' )
        
        
        number_seconds = len(all_seconds)
        for i, (second, num) in enumerate( zip(all_seconds, number_sources) ):
            if num < min_sources_in_second:
                print("skipping second", second, ". has too few sources")
                continue
            
            print(i, '/', number_seconds)
            delay_outputs, drift_outputs = CAL.process_second( second, num_runs, max_nfev=max_nfev )
            
            out_file.write(str(second) )
            out_file.write( '\n' )
            
            print()
            print()
            print('second', second,'done')
            for sname in CAL.antenna_names:
                
                ave_delay = 0.0
                err_delay = 0.0
                std_delay = 0.0
                N_delay = 0
                if sname in delay_outputs:
                    delays = delay_outputs[sname]
                    
                    ave_delay = np.average( delays )
                    N_delay = len( delays )
                    std_delay = np.std( delays )
                    err_delay = std_delay/np.sqrt(N_delay)
                
                ave_drift = 0.0
                err_drift = 0.0
                std_drift = 0.0
                N_drift = 0
                if sname in delay_outputs:
                    drifts = drift_outputs[sname]
                    
                    ave_drift = np.average( drifts )
                    N_drift = len( drifts )
                    std_drift = np.std( drifts )
                    err_drift = std_drift/np.sqrt(N_drift)
                

                print(sname, "N:", N_delay, ave_delay, '+/-', err_delay,  'STD:', std_delay)
                print("   ", ave_drift, '+/-', err_drift,  'STD:', std_drift)
                
                out_file.write( sname )
                out_file.write( ' ' )
                out_file.write( str(ave_delay) )
                out_file.write( ' ' )
                out_file.write( str(ave_drift) )
                out_file.write( ' ' )
                out_file.write( str(err_delay) )
                out_file.write( ' ' )
                out_file.write( str(err_drift) )
                out_file.write( '\n' )
            print()
            print()
    
    elif mode == 'test':
        ### TEST ONE SECOND ###
        second_to_test = None
        num_sources = 0
        for second, num in zip(all_seconds, number_sources):
            if num > num_sources:
                num_sources = num
                second_to_test = second
        
        
        delay_outputs, drift_outputs = CAL.process_second( second_to_test, num_runs, max_nfev=max_nfev )
            
        
        print()
        print()
        print('second', second,'done')
        for sname in CAL.antenna_names:
            
            ave_delay = 0.0
            err_delay = 0.0
            std_delay = 0.0
            N_delay = 0
            if sname in delay_outputs:
                delays = delay_outputs[sname]
                
                ave_delay = np.average( delays )
                N_delay = len( delays )
                std_delay = np.std( delays )
                err_delay = std_delay/np.sqrt(N_delay)
            
            ave_drift = 0.0
            err_drift = 0.0
            std_drift = 0.0
            N_drift = 0
            if sname in delay_outputs:
                drifts = drift_outputs[sname]
                
                ave_drift = np.average( drifts )
                N_drift = len( drifts )
                std_drift = np.std( drifts )
                err_drift = std_drift/np.sqrt(N_drift)
            

            print(sname, "N:", N_delay, ave_delay, '+/-', err_delay,  'STD:', std_delay)
            print("   ", ave_drift, '+/-', err_drift,  'STD:', std_drift)
    
    

