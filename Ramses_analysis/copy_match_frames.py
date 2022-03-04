import numpy as np
import os
import subprocess

Sink_49 = {'Projection_0': ['35900-36300'],
                'Projection_1': ['59300'],
                'Projection_2': ['58200-59500', '107400', '108000'],
                'Projection_3': ['35800', '58100-59600', '97700', '98700', '104300-104600'],
                'Projection_4': ['98600-98800', '104300'],
                'Projection_5': ['97800-97900', '98600'],
                'Projection_6': ['97100-97500'],
                'Projection_7': ['58000-58200', '104700-104800']}
Sink_91 = {'Projection_0': ['51800-52100'],
                'Projection_1': ['44800', '51900-52200', '55800-56000'],
                'Projection_2': ['50300-50800', '52000-52100'],
                'Projection_3': ['39500'],
                'Projection_4': ['49100'],
                'Projection_5': ['41400-41500'],
                'Projection_6': ['49100-49200', '50300'],
                'Projection_7': ['18800-48900']}
Sink_91HR = {'Projection_0': ['56800'],
                'Projection_1': ['53600'],
                'Projection_3': ['53000-53100'],
                'Projection_4': ['43500'],
                'Projection_5': ['45000', '52800', '53300'],
                'Projection_6': ['51600', '53100-53200', '56800'],
                'Projection_7': ['55200-55300']}
Sink_164 = {'Projection_1': ['45500'],
                'Projection_5': ['46800']}
                
dicts = [Sink_49, Sink_91, Sink_91HR, Sink_164]
dict_name = ['Sink_49', 'Sink_91', 'Sink_91HR', 'Sink_164']

dirs = ['/groups/astro/rlk/rlk/Movie_frames/Ramses/Sink_49/XY/1000AU/Obs_Comp/Image_Center_Primary/Time_Series/Centred_on_Secondary/Obs_threshold/', '/groups/astro/rlk/rlk/Movie_frames/Ramses/Sink_91/XY/1000AU/Obs_Comp/Image_Center_Primary/Time_Series/Obs_threshold/', '/groups/astro/rlk/rlk/Movie_frames/Ramses/Sink_91/High_resolution/XY/1000AU/Obs_Comp/Image_Center_Primary/Time_Series/Obs_threshold/', '/groups/astro/rlk/rlk/Movie_frames/Ramses/Sink_164/XY/1000AU/Obs_Comp/Image_Center_Primary/Time_Series/Low_resolution/Obs_threshold/']

dict_it = -1
for dict in dicts:
    dict_it = dict_it + 1
    for key in dict.keys():
        for match_time in dict[key]:
            if '-' in match_time:
                start_time = int(match_time.split('-')[0])
                end_time = int(match_time.split('-')[1])
                if end_time<start_time:
                    print('Error with programming match times!')
                    import pdb
                    pdb.set_trace()
                times = np.arange(start_time, end_time+100, 100)
            else:
                times = [int(match_time)]
            for time in times:
                call_string = 'rsync -vaz astro07-travel:'+dirs[dict_it] + 'movie_frame_' +("%06d" % (time/100))+'/p'+key[1:]+'.jpg '+str(dict_name[dict_it])+'/'+key+'/time_'+("%06d" % time)+'.jpg'
                call_string_rv = 'rsync -vaz astro07-travel:'+dirs[dict_it] + 'Zoom_in/movie_frame_' +("%06d" % (time/100))+'/p'+key[1:]+'_rv.jpg '+str(dict_name[dict_it])+'/'+key+'/time_'+("%06d" % time)+'_rv.jpg'
                
                subprocess.call(call_string, shell=True)
                subprocess.call(call_string_rv, shell=True)
