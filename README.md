radiation_operator.m    :    Main workhorse calculating fits of radiation operator
radiation_wrapper.m     :    Wrapper to pre- and post- process data of radiation_operator.m. Shows error for a range of weighting factors, so best can be selected.
readpltdata.m           :    Read data from plt_data/* (formatted ADAS data). The only necessary additional function if skipping radiation_wrapper.m.
Remaining *.m files     :    Called by either radiation_wrapper or radiation_operator
  create_rad_stucture.m :    Put radiation data into easy to use structure
  error_analysis.m      :    Check if fit matches desired criteria
  nuOfV.m               :    Nu in the operator
  rad_plot.m            :    plot data and fit
  vth.m                 :    v thermal
  scaled_mass.m         :    scaling factor for V0
  plt_data/*            :    Formatted ADAS data (human readable).
