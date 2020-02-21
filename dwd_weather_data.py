# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 10:30:04 2020

@author: Marian.Hummel
"""

class WeatherData:
    '''Create comprehensive weather data set with full hourly frequency.
    
    The goal is to achieve a comprehensive weather data set with the categories
    wind, sun, rain, clouds, temperature, humidity for as many common weather 
    stations as possible.
    
    Hourly weather data of the Deutscher Wetterdienst (DWD) are available per 
    category on the homepage as zip file. There are 2 files per category and
    measuring station (one with historical and one with current data).
    
    https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/
    
    Program can be started and controlled from the console.
    
        
    Methods
    ----------    
    list_zipfiles(list_of_urls, file_extension)
        Extracts all links of the required weather categories.
    
    filter_zipfiles()
        Filters out weather stations that do not exist for all desired weather
        categories and the historical and current data set.
        
    unpack_zipfiles(path_download_folder, exlucde_filter)
        Downloads the files temporarily and unpacks only the data set.
        
    import_weather_data(path)
        Loads all data records into the memory and names them automatically, 
        e.g. rain_hist_00161 for the historical values for rain of station
        00161. Each dataframe is saved as a single entry in a dictionary.
        
    create_city_dict()
        Create a dictionary with weather station ids and city names.
    
    merge_weather_data(data_to_merge)
        Vertically merge the historical and current dataframe per category and
        per station. These are then merged horizontally into one df per city.
        
    clean_weather_data(data_to_clean)
        Formats and renames columns. Creates a datetime object as index.
        
    interpolate_weather_data(false_values)
        Missing values labelled as -999 are interpolated.
        The dataframes consists of hourly data, but the frequency is not hourly 
        because there are missing observations on each day. These are also interpolated.
    '''    

    def __init__(self):
        
        pass
    
    
    ##############################################################
    # function to save urls of zip files
    def list_zipfiles(self, list_of_urls, file_extension='zip'):
        '''Extracts all links of the required weather categories.
        
        Parameters
        ----------
        list_of_urls : list
            List of links as str, eg.
            ['https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical']
        file_extension : str, optional
            Download only zip file, the default is 'zip'.

        Returns
        -------
        list
            List with urls of all zip files.
        '''

        # import
        from bs4 import BeautifulSoup
        import requests
        from tqdm import tqdm
        
        self.urls = list_of_urls
        
        self.zip_files_urls = list()
        for url in tqdm(self.urls, desc='save urls from server'):
                
            page = requests.get(url).text
            soup = BeautifulSoup(page, 'html.parser')
            self.zip_files_urls.append([url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(file_extension)])
        
        return self.zip_files_urls


    ##############################################################
    # Function to filter only the stations that are available for historical and current data in each weather category.
    def filter_zipfiles(self):
        '''Filters the full list of zip files.
        
        Returns
        -------
        list
            Only those weather stations for which all categories, historical 
            and current data are available.

        '''
     
        # all historical zip urls per weather feature
        zip_hist = [x for x in self.zip_files_urls if 'hist' in x[1]]
        # all historical zip urls in one list
        flattened_hist = [y for x in zip_hist for y in x] 
        
        # all current zip urls  per weather feature
        zip_curr = [x for x in self.zip_files_urls if 'akt' in x[1]]
        # all current zip urls in one list
        flattened_curr = [y for x in zip_curr for y in x] 

        # find common station ids per weather feature
        common_dict = dict()
        for index in range(len(zip_hist)):
            
            hist = list(map(lambda x: x.split('_')[-4], zip_hist[index]))
            curr = list(map(lambda x: x.split('_')[-2], zip_curr[index]))
            
            # common between hist and current for all weather features
            common_dict[index] = set.intersection(*map(set,[hist, curr]))
        
        # lowest common denominator of weather stations for all weather categories
        self.common = set.intersection(*map(set,[common_dict[0], common_dict[1], common_dict[2], common_dict[3], common_dict[4]]))
        
        hist_selected = [x for x in flattened_hist if x.split('_')[-4] in self.common]
        curr_selected = [x for x in flattened_curr if x.split('_')[-2] in self.common]
        
        # return list of filtered urls
        self.filtered_urls = list(hist_selected + curr_selected)
        
        return self.filtered_urls
                
    
    ##############################################################
    # function to download and extract zip files
    def unpack_zipfiles(self, path_download_folder, inclucde_filter = 'produkt_'):
        '''Downloads and unpacks the list or urls.
        
        By default only the table file is unpacked, all other files are ignored.

        Parameters
        ----------
        path_download_folder : str
            Path to target folder.
        inclucde_filter : str, optional
            Unpacks only files that start with inclucde_filter. The default is 'produkt_'.

        Returns
        -------
        None. Unzipped files are in destination folder.

        '''
       
        # import
        from io import BytesIO
        from urllib.request import urlopen
        from zipfile import ZipFile
        from tqdm import tqdm
 
        self.path = path_download_folder
        self.filter = inclucde_filter
        
        for url in tqdm(self.filtered_urls, desc='download and extract zipfiles'):
            with urlopen(url) as zipresp:
                with ZipFile(BytesIO(zipresp.read())) as zfile:
                    zfile.extractall(self.path ,
                                     members = (member for member in zfile.namelist() if member.startswith(self.filter))
                                     )

        print('\n downloaded and extracted data to: ', self.path)
        
        
    ##############################################################
    # function: load weather raw data
    def import_weather_data(self, path):
        """
        Loads raw weather data from as dataframes into a dictionary.
     
        
        Parameters
        ----------
        path : (relative) path to folder with raw weather data as string.
            E.g. './data/raw_weather/'
        
        Notes
        -----
        Two files exist for each weather feature and weather station. 
        One contains historical weather data up to 20XX and the other current data up to 20XX.
        The names of the data frames are formatted accordingly.
        
        For reasons of use and storage, only observations from the year 2015 
        onwards are retained. 
    
        
        Returns
        -------
        A dictionary with 2 data frames per station and weather feature.
        """
        import pandas as pd        
        import os
        from tqdm import tqdm
        
        self.path = path
        
        # dictionary with weather codes
        weather_codes = {'_tu' : 'air_temp',
                         '_ff' : 'wind',
                         '_rr' : 'rain',
                         '_sd' : 'sun',
                         '_n' : 'cloudiness'
                          }
            
        # create return dictionary
        self.weather_dict = dict()
        
        for file in tqdm(os.listdir(self.path), desc='load files from disc'):
        
            # create filename from the original name
            filename = os.fsdecode(file)
            filename = filename.replace('produkt','').replace('stunde_','').replace('.txt','').split('_')
            
            if filename[2].startswith('18') and filename[3].startswith('2018'):
                filename[2] = 'hist'
                
            elif filename[2].startswith('19') and filename[3].startswith('2018'):
                filename[2] = 'hist'
                
            elif filename[2].startswith('20') and filename[3].startswith('2018'):
                filename[2] = 'hist'
                
            elif filename[2].startswith('2018') and filename[3].startswith('2019'):
                filename[2] = 'current'
                
            elif filename[2].startswith('2018') and filename[3].startswith('2020'):
                filename[2] = 'current'
            
            else:
                print('Problem with', filename[4])
                
            filename.pop(3)
            filename = '_'.join(filename)
            
            # replace codes with weather feature
            for key in weather_codes.keys():
                filename = filename.replace(key, weather_codes[key])
                
            # load files and store filename in dictionary
            # drops unwanted columns directly in order to save memory!
            
            df = pd.read_csv(os.path.join(self.path, file), sep=';', usecols=[0,1,3,4])
            
            #dateparse = lambda x: pd.datetime.strptime(x, '%Y%m%d%H')
            #d = pd.read_csv(os.path.join(path, file), sep=';', usecols=[0,1,3,4], parse_dates=['MESS_DATUM'], date_parser = dateparse, index_col='MESS_DATUM', chunksize=50000)
            #data_frames[filename] = pd.concat([x[x.index.year >= 2015] for x in d])
             
            #df['date_time'] = pd.to_datetime(df['MESS_DATUM'], format='%Y%m%d%H')
            #df = df.set_index('date_time', sorted=True)
            #df = df[df.index.year >= 2015]
            
            # drop observations before 2015
            df = df.loc[df['MESS_DATUM'] > 2014123123]
            self.weather_dict[filename] = df
            
        return self.weather_dict
    
            
    ##############################################################    
    # function: create dictionary with station ids and city names 
    def create_city_dict(self):
        '''
        

        Returns
        -------
        None.

        '''
        
        # Certain attributes are required, therefore some methods must be reloaded if they have not been executed
        try:
            self.common
        except AttributeError:
            list_of_urls = ['https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/recent']
            self.list_zipfiles(list_of_urls)
            self.filter_zipfiles()
 
        # create dictionary containing all common station ids
        common_ids = list(self.common)
        f = open('./data/stations_id.txt','r')
        x = f.readlines()
        f.close()
        x = x[2:]
        
        self.city_codes = dict()
        for line in x:
            line = line.replace('  ',' ').replace('\n','').split()
            if line[0] in common_ids:
                station_id = line[0]
                station_name = line[6:]
                station_name.pop(-1)
                station_name = ' '.join(station_name)
                self.city_codes[station_id] = station_name       

                    
    ##############################################################
    # function: merge weather data
    def merge_weather_data(self, data_to_merge):
        """
        Merges all data frames into one.
        
        Parameters
        ----------    
        data : A dictionary containing all dataframes.
        
        Notes
        -----
        Two files exist for each weather feature and weather station. Therefore, 
        for each feature and station, the historical data must first be vertically 
        joined with the current data.
        
        These are then merged horizontally, duplicate entries are deleted and unused
        columns with the same name are dropped. The same applies to entries before 2015
        
        Returns
        -------
        One dataframe per city in dictionary
        """    
        import pandas as pd
        from functools import reduce
        from tqdm import tqdm
        
        # for passing already loaded data 
        self.weather_dict = data_to_merge
        
        # dictionary with one dataframe per city
        self.weather_merged = dict() 
        
        # list of problematic cities
        self.problem_cities = list()
        
        # Certain attributes are required, therefore some methods must be reloaded if they have not been executed
        try:
            self.common
        except AttributeError:
            list_of_urls = ['https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/recent']
            self.list_zipfiles(list_of_urls)
            self.filter_zipfiles()
        
        try:
            self.city_codes
        except AttributeError:
            self.create_city_dict()        
        
               
        # Create a data frame for each city that contains all weather categories. 
        for city in tqdm(self.city_codes, desc="merge dataframes into one per city"):
             
            # temporary dictionaries for vertical join
            dfs_wind = dict()
            dfs_sun = dict()
            dfs_cloudiness = dict()
            dfs_rain = dict()
            dfs_air_temp = dict()   
            
            # list of joined weather features per station
            weather_features = list()
            
            # final list of weather features to merge
            merge_list = list()
    
            # Check which items in dictionary have the same weather attributes 
            # concate historical and current datasets / vertical joining
            for (key, df) in self.weather_dict.items():
                
                # Check if keys are same, then add pair to new dictionary
                # if one pair per category (e.g. historical and current for rain)
                # are loaded, join
                
                if 'wind' in key and city in key:
                    dfs_wind[key] = df    
               
                    if len(dfs_wind) == 2:
                        wind = pd.concat(dfs_wind, sort = False)   
                        weather_features.append(wind)
            
                elif 'sun' in key and city in key:
                    dfs_sun[key] = df
                    
                    if len(dfs_sun) == 2:
                        sun = pd.concat(dfs_sun, sort = False)  
                        weather_features.append(sun)
            
                elif 'cloudiness' in key and city in key:
                    dfs_cloudiness[key] = df
                    
                    if len(dfs_cloudiness) == 2:
                        cloudiness = pd.concat(dfs_cloudiness, sort = False)   
                        weather_features.append(cloudiness)
            
                elif 'rain' in key and city in key:
                    dfs_rain[key] = df
                    
                    if len(dfs_rain) == 2:
                        rain = pd.concat(dfs_rain, sort = False)   
                        weather_features.append(rain)
            
                elif 'air_temp' in key and city in key:
                    dfs_air_temp[key] = df
                    
                    if len(dfs_air_temp) == 2:
                        air_temp = pd.concat(dfs_air_temp, sort = False)  
                        weather_features.append(air_temp)
            
                else:
                    continue
                
                      
            # pre processing for merge: clean overlapping periods, clean unwanted columns which are in all datasets
            # historical and current data set have overlapping periods
            
            for df in weather_features:
                df = df.droplevel(0)
                
                # drop duplicates
                df = df[~df['MESS_DATUM'].duplicated()]
                
                # one case eor doesn#t get dropped                
                df = df[df.columns.drop(list(df.filter(regex='eor')))]
                merge_list.append(df)
            
            try:
                df = reduce(lambda left,right: pd.merge(left,right,on=['MESS_DATUM', 'STATIONS_ID']), merge_list)
            
                # save final dataframe per city in dictionary
                self.weather_merged[city] = df
                
            except:
                print('Problem with:', city)
                
                self.problem_cities.append(city)
           
        return self.weather_merged
    

    ##############################################################
    # function: clean weather raw data  
    def clean_weather_data(self, data_to_clean):
        """
        Clean weather data and format column names
        
        Parameters
        ----------    
        data_to_clean : dict
            Dictionary with one dataframe per station.
        
        Notes
        -----
        Stations IDS will be converted to city names, weather feature codes will be renamed.
        Unused weather columns are deleted.
        
        Returns
        -------
        dict
            Cleaned dataframe per city in dictionary
        """
        
        import pandas as pd
        import tqdm
        
        # for passing already loaded data 
        self.weather_merged = data_to_clean
        
        # for saving cleaned data
        self.weather_cleaned = dict()
        
        # Certain attributes are required, therefore some methods must be reloaded if they have not been executed
        try:
            self.common
        except AttributeError:
            list_of_urls = ['https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/recent', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/historical', 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/recent']
            self.list_zipfiles(list_of_urls)
            self.filter_zipfiles()
        
        try:
            self.city_codes
        except AttributeError:
            self.create_city_dict()
        
        
        # iterate over dictionary with merged data sets and format data
        for (key,df) in tqdm.tqdm(self.weather_merged.items(), desc="format data"):
                        
            # make each element in the columns lowercase and remove whitespace
            df.columns = df.columns.str.lower().str.strip()
            
            # rename stations id according to city names
            df.rename(columns = {'stations_id' : 'stations_name'}, inplace = True)
            
            temp_city_codes = dict()
            for (k,v) in self.city_codes.items():
                k = int(k)
                temp_city_codes[k]=v
            
            df['stations_name'].replace(temp_city_codes, inplace = True)
            
            # rename columns of weather features with weather codes dictionary
            weather_codes = {'v_n' : 'cloudiness',
                             'rf_tu': 'rel_humidity',
                             'tt_tu' : 'air_temp',
                             'f': 'wind_speed',
                             'd': 'wind_direction',
                             'sd_so' : 'sunshine',
                             'r1': 'rain',
                            }
            
            df.rename(columns=weather_codes, inplace= True)   
        
            
            # create timestamp
            df['date_time'] = pd.to_datetime(df['mess_datum'], format='%Y%m%d%H')
            df = df.set_index('date_time')
            df.sort_index(ascending=True, inplace=True)
            try:
                df.drop(['mess_datum', 'v_n_i', 'rs_ind'], axis=1, inplace=True)
            except:
                pass
    
            self.weather_cleaned[key] = df                    
         
        return self.weather_cleaned
    
    
    ##############################################################
    # function: interpolate missing weather data
    
    def interpolate_weather_data(self, data_to_interpolate, false_values):
        '''Interpolate incorrect values and missing values due to incomplete 
        hourly frequency.
        
        Note
        ----
        Attention, this step should only be performed if interpolated data is wanted.
        

        Parameters
        ----------
        data: dict
            Dictionary with one dataframe per station.
        
        false_values : int
            Wrong value as integer, e.g. -999.
            

        Returns
        -------
        dict
            Interpolated data with full hourly frequency.
       
        '''
        import numpy as np
        import tqdm
        
        self.weather_cleaned = data_to_interpolate
        
        self.weather_final = dict()
        
        for (key,df) in tqdm.tqdm(self.weather_cleaned.items(), desc="interpolate data"):
        
            df = df.replace({false_values:np.NaN}, inplace=False)
            
            df = df.asfreq(freq='H')
            
            df = df.interpolate(method="time")
            
            df['stations_name'] = df['stations_name'].fillna(method="ffill")
            
            self.weather_final[key] = df
        
        return self.weather_final
    
    
# All information to initialize download, merging, cleaning and interpolation
def init_app():
    ### Part 1: download and extract raw data or load from savepoint
    
    # list of zip files to download (temp, clouds, rain, sun, wind)               
    list_of_urls = ['https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/recent',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/historical',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/cloudiness/recent' ,
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/historical',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/precipitation/recent',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/historical',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/sun/recent',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/historical',
                    'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/wind/recent']
        
    
    choice = input("The following start options are available for the further procedure: \
             \n(1) all links must be downloaded, unpacked and further processed (~11GB)\
             \n(2) only load from a savepoint 1 (~1.5GB, raw data in memory) \
             \n(3) load savepoint 2 (cleaned and interpolated! data ~0.5GB )?' \nPlease select (1/2/3/none) -> ")
    if choice == 1:
        # get list of urls
        DWD.list_zipfiles(list_of_urls)
        # filter list of urls to match only stations with full features and coverage
        DWD.filter_zipfiles()
        # download all zip files and extract to folder, keep only data from the year 2015+
        DWD.unpack_zipfiles('./data/raw_weather')
        # import each txt file to a dictionary of dataframes
        dfs_weather = DWD.import_weather_data('./data/raw_weather')
        # create savepoint1 ?
        if input("Create savepoint1 for later use (~1.5GB, imported raw data)? (y/n) -> ").lower().strip()[:1] == "y": 
            try:
                import hickle as hkl
                hkl.dump(dfs_weather, './data/weather/savepoint1.hkl')
                print('Savepoint 1 created')
            except ModuleNotFoundError:
                print('Not created, module not installed, see install.txt')
        # merge data
        dfs_weather_merged = DWD.merge_weather_data(dfs_weather)
        # clean and format dataframe per city
        dfs_weather_clean = DWD.clean_weather_data(dfs_weather_merged)
        # interpolate missing data per city
        dfs_weather_final = DWD.interpolate_weather_data(dfs_weather_clean, -999)
        print('All data merged, cleanded and interpolated. End.')
        # create final savepoint2 ?
        if input("Create savepoint2 for later use (~0.5GB, cleaned and interpolated data)? (y/n) -> ").lower().strip()[:1] == "y": 
            try:
                import hickle as hkl
                hkl.dump(dfs_weather_final, './data/weather/savepoint2.hkl')
                print('Savepoint 2 created')
            except ModuleNotFoundError:
                print('Not created, module not installed, see install.txt')
        # create one dataframe with all stations
        import pandas as pd
        weather = pd.concat(dfs_weather_final, axis=0)
        weather.to_csv('./data/weather.csv')
        print('final csv exported')
        
    elif choice == 2:
        try:        
            # module needs to be installed, see install.txt
            import hickle as hkl
            # savepoint path needs to be specified
            print('Please wait some time')
            dfs_weather = hkl.load('./data/weather/savepoint1.hkl')
            print('Savepoint loaded.')
            
        except OSError:
            print('dfs_weather.hkl does not exist. Please restart and download data.')   
        except ModuleNotFoundError:
            import sys
            sys.exit('Error: hickle module needs to be installed. Please see install.txt')
                
            
    elif choice == 3:
        try:        
            # module needs to be installed, see install.txt
            import hickle as hkl
            # savepoint path needs to be specified
            print('Please wait some time')
            dfs_weather_final = hkl.load('./data/weather/dfs_weather_final.hkl')
            print('Done.')
            
        except OSError:
            print('dfs_weather_final.hkl does not exist. Please restart and download data.')   
            
        except ModuleNotFoundError:
            import sys
            
            sys.exit('Error: hickle module needs to be installed, please see install.txt')
    else:
        print('Exit, need data to proceed. Please restart.')
        
        



###################################################       
if __name__ == "__main__":
    
    # create instance
    DWD = WeatherData()
    init_app()
    

