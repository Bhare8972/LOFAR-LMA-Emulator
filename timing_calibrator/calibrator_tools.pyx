#cython: language_level=3 
#cython: cdivision=True
#cython: boundscheck=True
#cython: linetrace=False
#cython: binding=False
#cython: profile=False


cimport numpy as np
import numpy as np
from libcpp cimport bool
from libc.math cimport sqrt, fabs, isfinite
from libc.stdlib cimport malloc, free

cdef double c_air = 299732445.9660432 
# cdef double c_air = 299792458.0/1.000293
cdef double c_air_inverse = 1.0/c_air

    
cdef class delay_fitter:
    
    cdef int num_station_delays 
    cdef int num_station_drifts 
    cdef int num_antennas 
    cdef int max_num_events
    cdef int current_num_events
    
    cdef double[:,:] antenna_locations # three nums per antenna
    cdef long[:] delay_indeces # one num per antenna. Index indicates which delay to use for antenna. negative indicates no extra delay
    cdef long[:] drift_indeces # one num per antenna. Index indicates which drift to use for antenna. negative indicates no drift
    cdef double[:] measurement_times # max_num_events*num_antennas
    cdef np.uint8_t[:] measurement_filter
    
    cdef int total_num_measurments
    
    cdef double drift_zero
    

    
    def __init__(self, np.ndarray[double, ndim=2] _antenna_locations, _max_num_events):
        self.num_antennas = _antenna_locations.shape[0]
        self.max_num_events = _max_num_events
        
        self.antenna_locations = _antenna_locations
        
        self.measurement_times = np.empty( _max_num_events*self.num_antennas, dtype=np.double )
        self.measurement_filter = np.zeros( _max_num_events*self.num_antennas, np.uint8 )
        
        self.current_num_events = 0
        self.total_num_measurments = 0
        self.drift_zero = 0.0
        
    def set_event(self, np.ndarray[double, ndim=1] arrival_times, np.ndarray[np.uint8_t, ndim=1] mask):
        filter_slice = self.measurement_filter[self.current_num_events*self.num_antennas : (self.current_num_events+1)*self.num_antennas]
        measurment_slice = self.measurement_times[self.current_num_events*self.num_antennas : (self.current_num_events+1)*self.num_antennas]

        
        cdef int ant_i
        for ant_i in range(self.num_antennas):
            if isfinite( arrival_times[ant_i] ):
                filter_slice[ant_i] = mask[ ant_i ]
                self.total_num_measurments += mask[ ant_i ]
            else:
                filter_slice[ant_i] = 0
                
            measurment_slice[ant_i] = arrival_times[ant_i] 
            
        self.current_num_events += 1
        
    def set_event_mask(self, event_i, station_i, np.uint8_t mask_value ):
        filter_slice = self.measurement_filter[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
        
        if filter_slice[station_i] and (not mask_value):
            self.total_num_measurments -= 1
        if (not filter_slice[station_i]) and mask_value:
            self.total_num_measurments += 1
            
        if mask_value:
            filter_slice[station_i] = 1
        else:
            filter_slice[station_i] = 0
        
        
    def set_delay_indeces(self, np.ndarray[long, ndim=1] _delay_indeces):
        self.num_station_delays = np.max( _delay_indeces ) + 1
        self.delay_indeces = _delay_indeces
        
    def set_drift_indeces(self, np.ndarray[long, ndim=1] _drift_indeces):
        self.num_station_drifts = np.max( _drift_indeces ) + 1
        self.drift_indeces = _drift_indeces
        
    def get_num_delays(self):
        return self.num_station_delays
        
    def get_num_drifts(self):
        return self.num_station_drifts
        
    def set_drift_zero(self, double drift_zero):
        self.drift_zero = drift_zero
        
    def reset(self):
        self.current_num_events = 0
        self.total_num_measurments = 0
        
    def firstOrder_objective_fun(self, np.ndarray[double , ndim=1] guess):
        cdef np.ndarray[double , ndim=1] ret = np.empty(self.total_num_measurments)
        
        cdef double[:] drifts = guess[ : self.num_station_drifts ]
        cdef double[:] delays = guess[ self.num_station_drifts : self.num_station_drifts+self.num_station_delays]
        cdef double[:] event_XYZTs = guess[ self.num_station_drifts+self.num_station_delays : ]
            
        cdef double X
        cdef double Y
        cdef double Z
        cdef double T
        
        
        cdef double dx
        cdef double dy
        cdef double dz
        cdef double dt
        
        cdef int event_i
        cdef int antenna_i
        cdef int delay_i
        cdef int drift_i
        cdef double total_delay
        cdef double antenna_drift
        cdef double measured_T
        cdef int output_i = 0
        cdef double[:] measurement_slice
        cdef np.uint8_t[:] filter_slice
        for event_i in range(self.current_num_events):
            measurement_slice = self.measurement_times[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
            filter_slice = self.measurement_filter[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
                        
            X = event_XYZTs[ event_i*4 + 0]
            Y = event_XYZTs[ event_i*4 + 1]
            Z = event_XYZTs[ event_i*4 + 2]
            Z = fabs(Z)
            T = event_XYZTs[ event_i*4 + 3]
            
            for antenna_i in range(self.num_antennas):
                if filter_slice[antenna_i]:
                    
                    dx = self.antenna_locations[antenna_i, 0] - X
                    dy = self.antenna_locations[antenna_i, 1] - Y
                    dz = self.antenna_locations[antenna_i, 2] - Z
                    
                    total_delay = 0.0
                    antenna_drift = 0.0
                    
                    delay_i = self.delay_indeces[ antenna_i ]
                    drift_i = self.drift_indeces[ antenna_i ]
                    if delay_i >= 0:
                        total_delay = delays[ delay_i ]
                        
                    if drift_i >= 0:
                        antenna_drift = drifts[ drift_i ]
                        
                    dt = sqrt(dx*dx + dy*dy + dz*dz)*c_air_inverse + T
                    
                    measured_T = measurement_slice[ antenna_i ]
                    dt -= measured_T - ( total_delay + antenna_drift*( measured_T-self.drift_zero ) )
                    
                    ret[output_i] = dt
                        
                    output_i += 1
                    
        return ret
        
        
    def full_objective_fun(self, np.ndarray[double , ndim=1] guess):
        cdef np.ndarray[double , ndim=1] ret = np.empty(self.total_num_measurments)
        
        cdef double[:] delays = guess[ : self.num_station_delays]
        cdef double[:] event_XYZTs = guess[ self.num_station_delays : ]
            
        cdef double X
        cdef double Y
        cdef double Z
        cdef double T
        
        cdef double dx
        cdef double dy
        cdef double dz
        cdef double dt
        
        cdef int event_i
        cdef int antenna_i
        cdef int delay_i
        cdef double total_delay
        cdef int output_i = 0
        cdef double[:] measurement_slice
        cdef np.uint8_t[:] filter_slice
        for event_i in range(self.current_num_events):
            measurement_slice = self.measurement_times[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
            filter_slice = self.measurement_filter[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
                        
            X = event_XYZTs[ event_i*4 + 0]
            Y = event_XYZTs[ event_i*4 + 1]
            Z = event_XYZTs[ event_i*4 + 2]
            Z = fabs(Z)
            T = event_XYZTs[ event_i*4 + 3]
            
            for antenna_i in range(self.num_antennas):
                if filter_slice[antenna_i]:
                    
                    dx = self.antenna_locations[antenna_i, 0] - X
                    dy = self.antenna_locations[antenna_i, 1] - Y
                    dz = self.antenna_locations[antenna_i, 2] - Z
                    
                    total_delay = 0.0
                    delay_i = self.delay_indeces[ antenna_i ]
                    if delay_i >= 0:
                        total_delay = guess[ delay_i ]
                        
                    dt = sqrt(dx*dx + dy*dy + dz*dz)*c_air_inverse + T + total_delay
                    dt -= measurement_slice[ antenna_i ]
                    
                    ret[output_i] = dt
                        
                    output_i += 1
                    
        return ret
    
    def sort_residuals(self, np.ndarray[double , ndim=1] residuals, out_dict, antenna_names):
        
        cdef int output_i = 0
        cdef int event_i
        cdef int antenna_i
        cdef double res
        cdef np.uint8_t[:] filter_slice
        for event_i in range(self.current_num_events):
            filter_slice = self.measurement_filter[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
            
            for antenna_i in range(self.num_antennas):
                if filter_slice[antenna_i]:
                    
                    res = residuals[ output_i ]
                    ant_name = antenna_names[ antenna_i ]
                    if ant_name not in out_dict:
                        out_dict[ ant_name ] = []
                    out_dict[ ant_name ].append( res )
                    
                    output_i += 1
                    
    def find_worse_fit(self, np.ndarray[double , ndim=1] residuals):
        
        cdef int event_i_out=0
        cdef int station_i_out=0
        cdef double worst_res=0
        
        cdef int event_i
        cdef int antenna_i
        cdef int output_i = 0
        cdef double res
        for event_i in range(self.current_num_events):
            filter_slice = self.measurement_filter[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
            
            for antenna_i in range(self.num_antennas):
                if filter_slice[antenna_i]:
                    res = residuals[ output_i ]
                    
                    if np.abs(res) > np.abs(worst_res):
                        worst_res = res
                        event_i_out = event_i
                        station_i_out = antenna_i
                    
                    output_i += 1
                    
        return event_i_out, station_i_out, worst_res
        
        
    def loc_objective_fun(self, np.ndarray[double , ndim=1] guess):
        cdef np.ndarray[double , ndim=1] ret = np.empty(self.total_num_measurments)
        
        cdef double[:] event_XYZTs = guess
            
        cdef double X
        cdef double Y
        cdef double Z
        cdef double T
        
        cdef double dx
        cdef double dy
        cdef double dz
        cdef double dt
        
        cdef int event_i
        cdef int antenna_i
        cdef int output_i = 0
        cdef double[:] measurement_slice
        cdef np.uint8_t[:] filter_slice
        for event_i in range(self.current_num_events):
            measurement_slice = self.measurement_times[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
            filter_slice = self.measurement_filter[event_i*self.num_antennas : (event_i+1)*self.num_antennas]
                        
            X = event_XYZTs[ event_i*4 + 0]
            Y = event_XYZTs[ event_i*4 + 1]
            Z = event_XYZTs[ event_i*4 + 2]
            Z = fabs(Z)
            T = event_XYZTs[ event_i*4 + 3]
            
            for antenna_i in range(self.num_antennas):
                if filter_slice[antenna_i]:
                    
                    dx = self.antenna_locations[antenna_i, 0] - X
                    dy = self.antenna_locations[antenna_i, 1] - Y
                    dz = self.antenna_locations[antenna_i, 2] - Z
                        
                    dt = sqrt(dx*dx + dy*dy + dz*dz)*c_air_inverse + T
                    dt -= measurement_slice[ antenna_i ]
                    
                    ret[output_i] = dt
                        
                    output_i += 1
                    
        return ret

                
        
    
    
    
    
    