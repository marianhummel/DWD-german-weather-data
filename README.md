# DWD-german-weather-data
Creates a comprehensive set of German weather data with full hourly frequency.
    
The aim is to obtain a comprehensive weather data set with the weather categories 
wind, sun, rain, clouds, temperature and humidity. All weather stations for which 
these measurements are available will be considered.

Note: The data source is the Deutscher Wetterdienst (DWD). Full information on
the data description can be found here
https://www.dwd.de/EN/ourservices/opendata/opendata.html and here
and in the folder of the respective observation category.

Hourly weather data of the Deutscher Wetterdienst (DWD) are available per 
category on the homepage as zip file. There are 2 files per category and
measuring station (one with historical and one with current data), see
https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/ .

Data is downloaded, unzipped, merged per weather station, cleaned and interpolated.

Attention, the last step was only carried out due to the need to obtain a full hourly 
frequency. The final data should therefore be analysed with appropriate caution.
