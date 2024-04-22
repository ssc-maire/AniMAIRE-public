import datetime as dt
import glob as gb
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import plotly.express as px


def convertSixStringIntoDatetime(TIME):
    hour = float(TIME[0:2])
    minute = float(TIME[2:4])
    second = float(TIME[4:6])

    return dt.timedelta(seconds = second + (60 * minute) + (60 * 60 * hour))

def getDatetimesForGLEDBdata(YYMMDD,TIME):

    date = dt.datetime.strptime(YYMMDD,"%y%m%d").replace(tzinfo=dt.timezone.utc)
    initialTime = convertSixStringIntoDatetime(TIME.split("-")[0])
    endTime = convertSixStringIntoDatetime(TIME.split("-")[1])

    initialDatetime = date + initialTime
    endDatetime = date + endTime    
    
    return (initialDatetime, endDatetime)

def importNMdataCSV(filePath):
    outputDF = pd.read_csv(filePath,header=None,skiprows=12,delimiter="\s+",skipfooter=20,dtype=str,engine="python").dropna()
    try:
        outputDF.columns = ["STATION","YYMMDD","SEC","TIME (UT)","CODE","UNCORR.","PRESS.","CORR.","% INC.","DETRENDED"]
    except:
        try:
            outputDF_copy = outputDF.copy()
            outputDF_copy[0] = outputDF_copy[0] + " " + outputDF_copy[1]
            outputDF_copy.drop(1,axis=1,inplace=True)
            outputDF_copy.columns = ["STATION","YYMMDD","SEC","TIME (UT)","CODE","UNCORR.","PRESS.","CORR.","% INC.","DETRENDED"]
            outputDF = outputDF_copy
        except:
            outputDF[9] = outputDF[8]
            outputDF[8] = outputDF[7]
            outputDF[7] = outputDF[6]

            outputDF[6] = outputDF[5].apply(lambda x:"-" + x.split("-")[1]) 
            outputDF[5] = outputDF[5].apply(lambda x:x.split("-")[0])
            
            outputDF.columns = ["STATION","YYMMDD","SEC","TIME (UT)","CODE","UNCORR.","PRESS.","CORR.","% INC.","DETRENDED"]
    
    outputDF["initialDatetime"] = outputDF.apply(lambda row:getDatetimesForGLEDBdata(row["YYMMDD"],row["TIME (UT)"])[0],axis=1)
    outputDF["endDatetime"] = outputDF.apply(lambda row:getDatetimesForGLEDBdata(row["YYMMDD"],row["TIME (UT)"])[1],axis=1)

    # getting file metadata
    with open(filePath,"r") as inputFile:
        fullFileContents = inputFile.readlines()
    
    outputDF["latitude"] = float(fullFileContents[0].split("LATITUDE")[1].split("LONGITUDE")[0])
    if ("c059" in filePath) or ("c060" in filePath) or ("c065" in filePath) or ("c071" in filePath):
        outputDF["longitude"] = float(fullFileContents[0].split("LONGITUDE")[1].split("ALTITUDE")[0])
    elif ("c070" in filePath) or ("c073" in filePath):
        raw_longitude = float(fullFileContents[0].split("LONGITUDE")[1].split("ALTITUDE")[0])
        if raw_longitude >= 0.0:
            output_longitude = raw_longitude
        elif raw_longitude < 0.0:
            output_longitude = raw_longitude + 360.0
        outputDF["longitude"] = output_longitude
    else:
        print("unknown GLE! Please check whether longitude is in longitude east format or not and add case to here.")
        raise Exception
    outputDF["altitude(m)"] = float(fullFileContents[0].split("ALTITUDE")[1].split("M")[0])
    outputDF["preincreaseCountRate(cts/s)"] = float(fullFileContents[4].split("PRE-INCREASE AVERAGE COUNTING RATE")[1].split("COUNTS PER SECOND")[0].replace(",",""))

    outputDF["CORR."] = outputDF["CORR."].apply(float)
    outputDF["UNCORR."] = outputDF["UNCORR."].apply(float)
    outputDF["% INC."] = outputDF["% INC."].apply(float)

    #print(float(fullFileContents[4].split("PRE-INCREASE AVERAGE COUNTING RATE")[1].split("COUNTS PER SECOND")[0].replace(",","")))
    
    return outputDF

def getAllStationDataForDatetime(inputDatetime):
    
    stationCountRateList = []
    stationFilePathList = gb.glob("investigationData/*.dat")
    if len(stationFilePathList) == 0:
        raise Exception("Error: could not find investigationData folder, or any files in investigationData folder!")
    
    for inputFilePath in stationFilePathList:
        fullImportedDF = importNMdataCSV(inputFilePath)
        data_at_timestamp = fullImportedDF[fullImportedDF["initialDatetime"] == inputDatetime]
        if not len(data_at_timestamp) == 0.0:
            stationCountRateList.append(data_at_timestamp)

    outputDF =  pd.concat(stationCountRateList).reset_index()
    outputDF["latitudeNearest5"] = 5 * ((outputDF["latitude"]/5).apply(int))
    outputDF["longitudeNearest5"] = 5 * ((outputDF["longitude"]/5).apply(int))
    outputDF["numerical increase"] = (outputDF["% INC."]/100) * outputDF["preincreaseCountRate(cts/s)"]

    return outputDF

def plotNMlocations():

    nmStationDF = getAllStationDataForDatetime(dt.datetime(year=2000,month=7,day=14,hour=10,minute=35,second=0))
    nmStationFig = px.scatter(nmStationDF,x="longitude",y="latitude",color="STATION")
    nmStationFig.update_layout(xaxis_range=[0,360])
    nmStationFig.update_layout(yaxis_range=[-90,90])
    nmStationFig.show()

    return nmStationFig


#############################################################################################################

def getDatetimesAndDFsForFilePattern(inputFilePattern, beginningEventDatetime, ignoreAltitude = None):

    listOfRelevantFiles = gb.glob(f"testingRuns/{inputFilePattern}")

    try:
        hoursAndMinutes = [[float(filePath.split("kft")[1].split("_")[0]),
                            float(filePath.split("_")[1].split(".csv")[0])] for filePath in listOfRelevantFiles]
    except IndexError:
        hoursAndMinutes = [[float(filePath.split("Alts")[1].split("_")[0]),
                            float(filePath.split("_")[1].split(".csv")[0])] for filePath in listOfRelevantFiles]

    relevantYear = beginningEventDatetime.year
    relevantMonth = beginningEventDatetime.month
    relevantDay = beginningEventDatetime.day

    listOfDatetimesToUse = [dt.datetime(year=relevantYear,
                                        month=relevantMonth,
                                        day=relevantDay,
                                        hour=int(time[0]),
                                        minute=int(time[1]),
                                        second=0)
                            for time in hoursAndMinutes]

    listOfDFs = []
    for index, inputFile in enumerate(listOfRelevantFiles):
        MishevSnapshot = pd.read_csv(inputFile)

        if ignoreAltitude is not None:
            MishevSnapshot = MishevSnapshot[~((MishevSnapshot["altitude (km)"] < ignoreAltitude + 0.1) & (MishevSnapshot["altitude (km)"] > ignoreAltitude - 0.1))]

        MishevSnapshot["Datetime"] = listOfDatetimesToUse[index]
        MishevSnapshot["MinutesSinceEventStart"] = MishevSnapshot["Datetime"].apply(lambda x:(x-beginningEventDatetime).seconds)
        listOfDFs.append(MishevSnapshot)
    return listOfDatetimesToUse,listOfDFs

def createSingleDoseMapPlot(inputFilePath,
                              beginningEventDatetime = dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                              selectAltitude = None):

    DFtoUse = pd.read_csv(inputFilePath)

    if selectAltitude is not None:
        DFtoUse = DFtoUse[round(DFtoUse["altitude (km)"],4) == selectAltitude]

    doseRateMap = px.scatter(DFtoUse, x="longitude",y="latitude",color="adose",
                            symbol_sequence=["square"],
                            range_y=[-90,90],
                            range_x=[0,360])
    
    doseRateMap.update_traces(marker={'size': 10})
    doseRateMap.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1, range=[-90,90]))
    doseRateMap.update_xaxes(range=[0,360])
    doseRateMap.update_yaxes(range=[-90,90])
    doseRateMap.update_layout(xaxis_scaleanchor="y")

    doseRateMap.update_layout(autosize=False,
                            width=800,
                            height=600)

    doseRateMap.show()

    return doseRateMap

def createAnimatedDoseMapPlot(inputFilePattern = "MishevFullTimeInterval12kft*.csv",
                              beginningEventDatetime = dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                              ignoreAltitude = None, selectAltitude = 11.28, dose_type="adose"):

    listOfDatetimesToUse, listOfDFs = getDatetimesAndDFsForFilePattern(inputFilePattern, beginningEventDatetime,ignoreAltitude=ignoreAltitude)

    interpFunc = interp1d([(x-beginningEventDatetime).seconds for x in listOfDatetimesToUse],
            [listOfDFs[index][dose_type] for index in range(0,len(listOfDFs))],
            axis=0,kind='linear')

    max_minutes_to_use = (np.sort(listOfDatetimesToUse)[-1] - np.sort(listOfDatetimesToUse)[0]).seconds / 60

    interpedListOfDFs = []
    templateDF = listOfDFs[0]
    for minutes in np.linspace(0,max_minutes_to_use,int(max_minutes_to_use/5)+1):
        newInterpDF = templateDF.copy()
        newInterpDF[dose_type] = interpFunc(minutes*60)
        newInterpDF["MinutesSinceEventStart"] = minutes
        interpedListOfDFs.append(newInterpDF)
    
    #print(len(interpedListOfDFs))

    fullConcattedDF = pd.concat(interpedListOfDFs).sort_values(by=['Datetime','latitude','longitude'])
    fullConcattedDF = fullConcattedDF[round(fullConcattedDF["altitude (km)"],4) == selectAltitude]

    # getting NM data for plotting ####################################
    nmStationDF = pd.concat([getAllStationDataForDatetime(NMdatetime) for NMdatetime in listOfDatetimesToUse])
    nmStationDF["MinutesSinceEventStart"] = (nmStationDF["initialDatetime"] - listOfDatetimesToUse[0]).apply(lambda x:x.seconds / 60.0)
    nmStationDF.rename(columns={"longitude":"NMlongitude","latitude":"NMlatitude"})

    #NMscatter = px.scatter(nmStationDF,x="longitude",y="latitude",color="% INC.",animation_frame="MinutesSinceEventStart")
    ###################################################################

    mergedNMandMAIREdata = pd.merge(fullConcattedDF,nmStationDF,how="inner",on="MinutesSinceEventStart")

    doseRateMap = px.scatter(fullConcattedDF, x="longitude",y="latitude",color=dose_type,
                            symbol_sequence=["square"],
                            animation_frame="MinutesSinceEventStart",
                            range_y=[-90,90],
                            range_x=[0,360])

    # doseRateMap = px.scatter(mergedNMandMAIREdata, x=["longitude","NMlongitude"],y=["latitude","NMlatitude"],color=["adose","% INC."],
    #                         symbol_sequence=["square"],
    #                         animation_frame="MinutesSinceEventStart",
    #                         range_y=[-90,90],
    #                         range_x=[0,360])

    doseRateMap.update_traces(marker={'size': 10})
    doseRateMap.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1, range=[-90,90]))
    doseRateMap.update_xaxes(range=[0,360])
    doseRateMap.update_yaxes(range=[-90,90])
    doseRateMap.update_layout(xaxis_scaleanchor="y")

    doseRateMap.update_layout(autosize=False,
                            width=800,
                            height=600)
    #doseRateMap.update_layout(xaxis_range=[0,360],yaxis_range=[-90,90])

    #combinedOutputPlot = go.Figure(data=doseRateMap.data + NMscatter.data,layout=doseRateMap.layout)

    

    #doseRateMap.show()
    #combinedOutputPlot.show()

    #return (doseRateMap, listOfDatetimesToUse)
    return doseRateMap


###################################################################

def getCountRatesAndDoses(simulatedGCRDF,simulatedGLEDF,relevantTimestamp):

    nmCountRatesToUse = getAllStationDataForDatetime(relevantTimestamp)
    nmCountRatesToUse = nmCountRatesToUse[nmCountRatesToUse["altitude(m)"] < 300]
    nmCountRatesToUse = nmCountRatesToUse[nmCountRatesToUse["STATION"] != "INUVIK"]
    nmCountRatesToUse.reset_index(inplace=True,drop=True)

    outputDoseSeriesList = []
    for index, row in nmCountRatesToUse.iterrows():
        latitude = row["latitudeNearest5"]
        longitude = row["longitudeNearest5"]
        altitude = simulatedGLEDF["altitude (km)"].iloc[(simulatedGLEDF["altitude (km)"] - (row["altitude(m)"]/1000)).abs().idxmin()]
        #print(altitude)
        #print(simulatedGLEDF.head())

        GLEsimulatedDoseRateSeries = simulatedGLEDF[
                                            (simulatedGLEDF["latitude"] == latitude) &
                                            (simulatedGLEDF["longitude"] == longitude) &
                                            (simulatedGLEDF["altitude (km)"] == altitude)]
        GCRsimulatedDoseRateSeries = simulatedGCRDF[
                                            (simulatedGCRDF["latitude"] == latitude) &
                                            (simulatedGCRDF["longitude"] == longitude) &
                                            (simulatedGCRDF["altitude (km)"] == altitude)]     
        if len(GCRsimulatedDoseRateSeries) == 0:                 
            GCRsimulatedDoseRateSeries = simulatedGCRDF[
                                                (simulatedGCRDF["latitude"] == latitude) &
                                                (simulatedGCRDF["longitude"] == longitude) &
                                                (simulatedGCRDF["altitude (km)"] > altitude - 0.1) & 
                                                (simulatedGCRDF["altitude (km)"] < altitude + 0.1)]

        GCRsimulatedDoseRateSeries.reset_index(inplace=True)
        GLEsimulatedDoseRateSeries.reset_index(inplace=True)

        #print("+++++++++++++++++++++++++++++++++++++")
        #print(simulatedGCRDF[(simulatedGCRDF["altitude (km)"] == altitude)])
        #print(GCRsimulatedDoseRateSeries[["adose","tn1","tn2","tn3"]])
        #print(GLEsimulatedDoseRateSeries[["adose","tn1","tn2","tn3"]])
        #print(GLEsimulatedDoseRateSeries/GCRsimulatedDoseRateSeries)
        #print("+++++++++++++++++++++++++++++++++++++")

        simulatedDoseRateSeries = 100 * ((GLEsimulatedDoseRateSeries/GCRsimulatedDoseRateSeries))
        if len(simulatedDoseRateSeries) == 0:
            print(simulatedDoseRateSeries)
            print(latitude,longitude)

        #print(simulatedDoseRateSeries)
        
        outputDoseSeriesList.append(simulatedDoseRateSeries[["adose","tn1","tn2","tn3"]])
    fullOutputDoseDF = pd.concat(outputDoseSeriesList).reset_index()
    #print(len(fullOutputDoseDF))
    #print(len(nmCountRatesToUse))

    countRatesAndDosesGeneral = pd.concat([nmCountRatesToUse,fullOutputDoseDF],axis=1)
    countRatesAndDosesGeneral["logAltitudeInm"] = np.log(countRatesAndDosesGeneral["altitude(m)"])
    countRatesAndDosesGeneral.replace([np.inf, -np.inf], np.nan, inplace=True)
    return countRatesAndDosesGeneral

def getNMdataForAllDatetimes(inputGCRDF, inputFilePattern = "MishevFullTimeInterval12kft*.csv", 
                             beginningEventDatetime = dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                             ignoreAltitude = None):
    listOfDatetimesToUse, listOfGLEDFs = getDatetimesAndDFsForFilePattern(inputFilePattern = inputFilePattern,
                                beginningEventDatetime = beginningEventDatetime, ignoreAltitude=ignoreAltitude)
    #print(listOfDatetimesToUse)

    listOfNMscatterData = get_neutron_monitor_data(inputGCRDF, beginningEventDatetime, listOfDatetimesToUse, listOfGLEDFs)
    
    fullNMscatterData = pd.concat(listOfNMscatterData).sort_values(by=['MinutesSinceEventStart','STATION'])

    return fullNMscatterData

def get_neutron_monitor_data(inputGCRDF, beginningEventDatetime, listOfDatetimesToUse, listOfGLEDFs):
    listOfNMscatterData = []
    for index, inputGLEDF in enumerate(listOfGLEDFs):
        specificDatetime = listOfDatetimesToUse[index]
        newNMscatterData = getCountRatesAndDoses(inputGCRDF,
                        inputGLEDF,
                        specificDatetime)
        newNMscatterData["MinutesSinceEventStart"] = (specificDatetime - beginningEventDatetime).seconds / 60
        #print(specificDatetime)
        listOfNMscatterData.append(newNMscatterData)
    return listOfNMscatterData

def plotNMdataForAllDatetimes(inputGCRDF, inputFilePattern = "MishevFullTimeInterval12kft*.csv", 
                             beginningEventDatetime = dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                             ignoreAltitude = None,**xargs):
    
    fullNMscatterData = getNMdataForAllDatetimes(inputGCRDF, inputFilePattern = inputFilePattern, 
                             beginningEventDatetime = beginningEventDatetime,
                             ignoreAltitude = ignoreAltitude)

    #print(fullNMscatterData)

    fullScatterPlot = px.scatter(fullNMscatterData,x="% INC.",y="tn3",
                                color="STATION",
                                animation_frame="MinutesSinceEventStart",
                                #range_y=[-0.005,0.015],
                                #range_y=[-50,350],
                                #range_x=[-5,60],
                                width=800,
                                height=800,
                                **xargs)

    fullScatterPlot.add_scatter(x=np.linspace(-5,60,200),y=np.linspace(-5,60,200))

    # fullScatterPlot.update_yaxes(
    #                             scaleanchor="x",
    #                             scaleratio=1,
    #                             )

    return fullScatterPlot

def plotNMdataForAllDatetimesSquashed(inputGCRDF, inputFilePattern = "MishevFullTimeInterval12kft*.csv", 
                             beginningEventDatetime = dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                             ignoreAltitude = None, **xargs):
    
    fullNMscatterData = getNMdataForAllDatetimes(inputGCRDF, inputFilePattern = inputFilePattern, 
                             beginningEventDatetime = beginningEventDatetime,
                             ignoreAltitude = ignoreAltitude)

    fullScatterPlot = px.scatter(fullNMscatterData,x="% INC.",y="tn3",
                                color="STATION",
                                #animation_frame="MinutesSinceEventStart",
                                #color_continuous_scale="picnic",
                                width=800,
                                height=800,
                                **xargs)
    # fullScatterPlot.update_yaxes(
    #                             scaleanchor="x",
    #                             scaleratio=1,
    #                             )
    fullScatterPlot.add_scatter(x=np.linspace(-5,60,200),y=np.linspace(-5,60,200))

    fullScatterPlot.show()

def plotNMdataForSingleDatetime(inputGCRDF, inputFilePattern = "MishevFullTimeInterval12kft*.csv", 
                             beginningEventDatetime = dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                             MinutesSinceEventStart=0.0,
                             ignoreAltitude = None, **xargs):
    
    fullNMscatterData = getNMdataForAllDatetimes(inputGCRDF, inputFilePattern = inputFilePattern, 
                             beginningEventDatetime = beginningEventDatetime,
                             ignoreAltitude = ignoreAltitude)

    fullNMscatterData = fullNMscatterData.query(f"MinutesSinceEventStart == {MinutesSinceEventStart}")

    fullScatterPlot = px.scatter(fullNMscatterData,x="% INC.",y="tn3",
                                color="STATION",
                                #animation_frame="MinutesSinceEventStart",
                                #color_continuous_scale="picnic",
                                width=800,
                                height=800,
                                **xargs)
    # fullScatterPlot.update_yaxes(
    #                             scaleanchor="x",
    #                             scaleratio=1,
    #                             )
    fullScatterPlot.add_scatter(x=np.linspace(-5,60,200),y=np.linspace(-5,60,200))

    #fullScatterPlot.show()

    return fullScatterPlot