# LOFAR-LMA-Emulator
Code needed to emulate the function of a LMA using LOFAR data, and to calibrate the timing offsets of LMA data

Technically requires LoLIM to be installed as well ( https://github.com/Bhare8972/LOFAR-LIM )

LMA EMULATOR:

This highly relies upon LOLIM.

folder LMA_emulator requires essential files for emulating the operation of a LMA, including the different windowing techniqus
files:
  write_lma_data_10us.c
    this converts the output of the text output of the LMA emulator into 1 us binary format for procssing by the LMA algorithm
  
  run_traditionalWindowing_emulator_example.py
    an example of how to use the LMA emulator. Note that modifications are probably necisary to function. This example uses the traditional 80 us windowing.
    
  LMA_window_data.py
    the LMA emulator code for traditional windowing of any length. Some code is imported and re-used by other windowing codes.
    
  LMA_runThreshWindow_data.py
    same for floating threshold
    
  LMA_binnedWindow_data.py
    natural threshold (the name changed a few times)
    
  LMA_slidingWindow_data.py
    Non-Alligned Window
    
    
  steps to running the LMA-emulator:
    1) window the LOFAR data (use LMA_window_data.py,  LMA_runThreshWindow_data.py,  LMA_binnedWindow_data.py,  or  LMA_slidingWindow_data.py for the appropriate windowing technique)
    2) use write_lma_data_10us.c to convert the text output into LMA binary data
    3) process the binary with LMA 10 us processing
    
 
TIMING CALIBRATION
This technically only needs "read_LMA" from LoLIM. So it would to be possible to run this without fully installing LoLIM.

folder timing_calibrator contains necisary code to calibration the timing of LMA data.
Files:
  calibrator_tools.pyx
    cython code needed to run the calibrator
  setup_calibrator_tools.py
    python code to build the cython code
  build_calibrator_tools.sh
    one bash line to run the python code that builds the cython code. Nothing beats makeing a program easier to use than glueing four different langueges together.
    
  LMA_calibrator.py
    the actual calibrator. This is the file to run. Has lots of nice settings in the code itself.
    
  LMA_calibrator_FirstOrder.py
    very similar. But also fits first order of clock drift over each second. 
    
  


