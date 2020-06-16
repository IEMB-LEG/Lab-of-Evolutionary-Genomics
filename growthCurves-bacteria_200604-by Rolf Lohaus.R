# version 1.3.1h-200604 Scripted by Rolf Lohaus
# version 1.3g-120411


library(lattice)
library(reshape)


# convert time data from ddd:hh:mm:ss to seconds
seconds <- function(time) {
    t = as.numeric(noquote(unlist(strsplit(sub("\\.", ":", paste(time)), split=':', fixed=T))))
    sec = NULL
    if (length(t) == 4) {
        sec = t[1]*24*60*60+t[2]*60*60+t[3]*60+t[4]
    } else {
        if (length(t) == 3) {
            sec = t[1]*60*60+t[2]*60+t[3]
        } else {
            if (length(t) == 2) {
                sec = t[1]*60+t[2]
            } else {
                stop(paste0("Unsupported time format found in [", time, "], needs to conform to ddd:hh:mm:ss or subset thereof"))
            }
        }
    }

    return(sec)
}

minutesAndSeconds <- function(seconds) {
    return(c(seconds%/%60, seconds%%60))
}

avgSmooth <- function(opticalDensities, windowSize, log=T) {
    if (length(opticalDensities) <= 1) {
        return(0)
    }

    if (log) {
        opticalDensities = log(opticalDensities)
    }
    offset = (windowSize-1)/2
    expandedODs = c(rep(opticalDensities[1], offset), opticalDensities, rep(opticalDensities[length(opticalDensities)], offset))

    rates = NULL
    for (t in (1:length(opticalDensities))+offset) {
        mean = mean(expandedODs[(t-offset):(t+offset)])
        rates = c(rates, mean)
        # message(coef(m)[2])
    }

    return(rates)
}

# calculate derivative of parameter 'opticalDensities' using symmetrical window around data points
ratesInterval <- function(opticalDensities, times, windowSize, log=T) {
    if (length(opticalDensities) <= 1) {
        return(0)
    }

    if (log) {
        opticalDensities = log(opticalDensities)
    }
    offset = (windowSize-1)/2
    expandedODs = c(rep(opticalDensities[1], offset), opticalDensities, rep(opticalDensities[length(opticalDensities)], offset))

    deltaTime = times[2]-times[1]
    expandedTimes = c(seq(times[1]-deltaTime, times[1]-offset*deltaTime, -deltaTime), times, seq(times[length(times)]+deltaTime, times[length(times)]+offset*deltaTime, deltaTime))

    rates = NULL
    for (t in (1:length(opticalDensities))+offset) {
        m = lm(expandedODs[(t-offset):(t+offset)] ~ expandedTimes[(t-offset):(t+offset)])
        rates = c(rates, coef(m)[2])
        # message(coef(m)[2])
    }

    return(rates)
}

# calculate derivative of parameter 'opticalDensities' using 'windowLeft' and 'windowRight' points around data points
ratesOffsets <- function(opticalDensities, times, windowLeft, windowRight, log=T) {
    if (length(opticalDensities) <= 1) {
        return(0)
    }

    if (log) {
        opticalDensities = log(opticalDensities)
    }
    # message(offset)
    expandedODs = c(rep(opticalDensities[1],windowLeft), opticalDensities, rep(opticalDensities[length(opticalDensities)],windowRight))
    # print(expandedODs)

    deltaTime = times[2]-times[1]
    expandedTimes = c(seq(times[1]-deltaTime,times[1]-windowLeft*deltaTime,-deltaTime), times, seq(times[length(times)]+deltaTime,times[length(times)]+windowRight*deltaTime,deltaTime))

    slopes = NULL
    for (i in (1:length(opticalDensities))+windowLeft) {
        # message(i)
        m = lm(expandedODs[(i-windowLeft):(i+windowRight)] ~ expandedTimes[(i-windowLeft):(i+windowRight)])
        slopes = c(slopes, coef(m)[2])
    }

    return(slopes)
}


normalize <- function(x) {
    return((x-min(x))/(max(x)-min(x)) + (min(x)/(max(x)-min(x))))
}


isTwoResourceWell <- function(oneResourceWells, twoResourceWells, well) {
    if (is.null(twoResourceWells) || is.na(twoResourceWells) || length(twoResourceWells) == 0) {
        return(FALSE)
    }
    if (toupper(twoResourceWells[1]) == "ALL") {
        return(TRUE)
    }
    if (toupper(twoResourceWells[1]) == "NONE") {
        return(FALSE)
    }
    if (toupper(oneResourceWells[1]) == "ALL") {
        return(TRUE)
    }

    if (well %in% twoResourceWells) {
        return(TRUE)
    }

    rowCol = substring(well,1:2,c(1,3))
    if (rowCol[1] %in% twoResourceWells) {
        return(TRUE)
    }
    if (rowCol[2] %in% twoResourceWells) {
        return(TRUE)
    }

    if (length(oneResourceWells) != 0) {
        if (well %in% oneResourceWells) {
            return(FALSE)
        }

        rowCol = substring(well,1:2,c(1,3))
        if (rowCol[1] %in% oneResourceWells) {
            return(FALSE)
        }
        if (rowCol[2] %in% oneResourceWells) {
            return(FALSE)
        }

        return(TRUE)
    }

    return(FALSE)
}




rankPeaks <- function(od, rate, rate1, rate2) {
    return(1.5*od + (8*rate)^3 + 32*(1-abs(rate1)) - 4*rate2)
}


rankValley <- function(time, od, rate, rate1, rate2) {
    weightTime = 2
    weightRate = 2
    weightRate1 = 10
    weightRate2 = 7

    od1 = min(od) + (max(od) - min(od)) * 0.15
    od2 = min(od) + (max(od) - min(od)) * 0.75
    time1Index = which(od>od1)[1]
    time4Index = which(od>od2)[1]
    time1 = time[time1Index]
    time4 = time[time4Index]
    time2 = time1+(time4-time1)*0.2
    time3 = time1+(time4-time1)*0.5

    return(weightTime*sapply(time, rankValleyTime, time1, time2, time3, time4, max(time)) + weightRate*sapply(rate, rankValleyRate) + weightRate1*(1-abs(rate1)) + weightRate2*rate2)
}

rankValleyTime <- function(t, time1, time2, time3, time4, maxT) {
    if (t < time1 || t >= time4) {
        return(-0.4 - 1.2*t/maxT)
    }
    if (t < time2) {
        return((t-time1)/(time2-time1)-0.4 - 1.2*t/maxT)
    }
    if (t < time3) {
        return(0.6 - 1.2*t/maxT)
    }
    if (t < time4) {
        return(0.6-(t-time3)/(time4-time3) - 1.2*t/maxT)
    }
}

rankValleyRate <- function(x) {
    if (x >= 0) {
        return(1-(abs(x))^2)
    } else {
        return(1-6*(abs(x))^2)
    }
}





####### following 2 functions specific for transforming "08202009-2.txt" data set from Selwyn #######

threeCharWell <- function(wellName) {
    if (nchar(wellName) == 2) {
        return(paste(substr(wellName,1,1), substr(wellName,2,2), sep="0"))
    } else {
        return(wellName)
    }
}
createMissingWell <- function(well, times, od) {
    tmp = NULL
    for (t in times) {
        tmp = rbind(tmp, data.frame(Time=t, OD=od, Well=well))
    }
    return(tmp)
}






####### main analyses and plotting code starts here #######

data = NULL
newdata = NULL
wellNames = NULL

# baseODs <- NULL
# lags <- NULL
# ranksValley <- NULL
# ranksPeak1 <- NULL
# ranksPeak2 <- NULL
# valleys <- NULL
# peak1s <- NULL
# peak2s <- NULL
# expGrowths <- NULL
# models <- list()



analyzeGrowthCurves <- function(dataFileName, singleWell=NULL, oneResourceWells=c(), twoResourceWells=c("ALL"), plotDetail=F, reload=T, smoothWindowSize=7, maxOD=0.45, ODTicks=c(0.1,0.2,0.4), maxTimeHours=15, timeTicks=c(0,5,10,15), rSquaredTreshold=0.85, maxWells=96, selwyn=F) {
    if (!is.null(twoResourceWells) && !is.na(twoResourceWells) && length(twoResourceWells) != 0 &&
        !is.null(oneResourceWells) && !is.na(oneResourceWells) && length(oneResourceWells) != 0)
    {
        if (toupper(oneResourceWells[1]) == "ALL" && toupper(twoResourceWells[1]) == "ALL") {
            stop("ERROR: parameters oneResourceWells and twoResourceWells can't both be set to \"ALL\"")
        }
        if (toupper(oneResourceWells[1]) == "NONE" && toupper(twoResourceWells[1]) == "NONE") {
            stop("ERROR: parameters oneResourceWells and twoResourceWells can't both be set to \"NONE\"")
        }
    }

    ptm <- proc.time()

    # constants
    modelCol = rgb(0.2,0.9,0.2)
    modelColLight = rgb(0.6,1.0,0.6)
    maxTimeSeconds = maxTimeHours*60*60
    maxX = maxTimeSeconds*1.075
    minX = 0-maxTimeSeconds*0.075
    xTicks = timeTicks*60*60
    if (selwyn) {
        minY = 0.025
        maxY = 0.8
        yTicks = c(0.03, 0.3, 0.7)
        minOD = 0.075
        wellNameY = 0.55
        twoResourceY = 0.35
        lagValueY = 0.052
        lagValue1Y = 0.078
        lagValue2Y = 0.055
        slopeValueY = 0.035
        slopeValue1Y = 0.052
        slopeValue2Y = 0.035
        pY = 0.63
        residualsY = 0.45

        ratesYLim = c(-0.004,0.01)
    } else {
        minOD = 0.125
        yTicks = ODTicks
        minY = 0.08
        maxY = maxOD * 1.1
        logMinY = log(minY)
        logMaxY = log(maxY)
        ySize = logMaxY - logMinY
        wellNameY = exp(logMinY + 0.9 * ySize)
        twoResourceY = exp(logMinY + 0.78 * ySize)
        lagValue1Y = exp(logMinY + 0.46 * ySize)
        lagValue2Y = exp(logMinY + 0.35 * ySize)
        slopeValue1Y = exp(logMinY + 0.21 * ySize)
        slopeValue2Y = exp(logMinY + 0.1 * ySize)
        lagValueY = exp(logMinY + 0.21 * ySize)
        slopeValueY = exp(logMinY + 0.1 * ySize)
        pY = 0.63
        residualsY = 0.48

        ratesYLim = c(-0.55,0.95)
        ratesYTicks = c(-0.4, 0, 0.4, 0.8)
    }



    dataFile = unlist(strsplit(dataFileName, split=".txt", fixed=T))

    wellsWithZeros = NULL
    if (reload) {
        baseODs <<- NULL
        lags <<- NULL
        ranksValley <<- NULL
        ranksPeak1 <<- NULL
        ranksPeak2 <<- NULL
        valleys <<- NULL
        peak1s <<- NULL
        peak2s <<- NULL
        expGrowths <<- NULL
        models <<- list()
        problems <<- NULL

        message("Loading data from [", dataFileName, "]")
        if (!selwyn) {
            # data <<- read.table(dataFileName, header=T, skip=2, fill=T, fileEncoding = "latin1")
            data <<- read.csv(dataFileName, sep="\t", skip=2, header=T, fill=T, na.strings=c("", " ", "NA"), fileEncoding = "latin1")
            # message("\nHead:")
            # print(head(data, n=3))
            # message("\nTail:")
            # print(tail(data, n=5))

            # try to remove empty rows at end of file
            emptyRows=which(is.na(data$H12))
            if (length(emptyRows) != 0 && length(emptyRows) != nrow(data)) {
                # message("\nFound ", length(emptyRows), " empty rows in well H12")
                # print(emptyRows)
                # data <<- read.table(dataFileName, header=T, skip=2, fill=T, nrows=emptyRows[1]-1, fileEncoding = "latin1")
                data <<- read.csv(dataFileName, sep="\t", skip=2, header=T, fill=T, nrows=emptyRows[1]-1, na.strings=c("", " ", "NA"), fileEncoding = "latin1")
            }

            names(data)[1] <<- "Time"
            names(data)[2] <<- "Temperature"

            # try to remove empty rows at end of file
            emptyRows=which(is.na(data$Temperature))
            if (length(emptyRows) != 0 && length(emptyRows) != nrow(data)) {
                # message("\nFound ", length(emptyRows), " empty rows in column Temperature")
                # print(emptyRows)
                data <<- data[1:emptyRows[1]-1,]
            } else {
                emptyRows=which(is.na(data$B1))
                if (length(emptyRows) != 0 && length(emptyRows) != nrow(data)) {
                    # message("\nFound ", length(emptyRows), " empty rows in well B1")
                    # print(emptyRows)
                    data <<- data[1:emptyRows[1]-1,]
                } else {
                    emptyRows=which(is.na(data$H12))
                    if (length(emptyRows) != 0 && length(emptyRows) != nrow(data)) {
                        # message("\nFound ", length(emptyRows), " empty rows in well H12")
                        # print(emptyRows)
                        data <<- data[1:emptyRows[1]-1,]
                    }
                }
            }
            # message("\nTail:")
            # print(tail(data, n=5))

			# convert time data to seconds
            data$Time <<- sapply(data$Time, seconds)
            # message("\nTail:")
            # print(tail(data, n=5))


            message("\nPre-processing wells")

            message("\nSmoothing OD values...")

            if (is.null(singleWell)) {
                firstWellName = "A1"
            } else {
                firstWellName = singleWell
            }

            if (length(which(data[[firstWellName]]==0)) > 0) {
                data[[firstWellName]][which(data[[firstWellName]]==0)] = 0.001
                wellsWithZeros = c(wellsWithZeros, firstWellName)
                message(sprintf("***  %-3s contains zero-value data points  ***", firstWellName))
            }
            if (length(which(is.na(data[[firstWellName]]))) == nrow(data)) {
                data[[firstWellName]][which(is.na(data[[firstWellName]]))] = 0.0001
                wellsWithZeros = c(wellsWithZeros, firstWellName)
                message(sprintf("***  %-3s contains no data                 ***", firstWellName))
            }
            smooth = data.frame(Time=data$Time, OD=data[[firstWellName]], SmoothedOD=avgSmooth(data[[firstWellName]], smoothWindowSize,log=F), Well=firstWellName)

            if (is.null(singleWell)) {
                for (w in 2:(min(maxWells, length(data)-2))) {
                    wellName = names(data)[w+2]
                    # message(wellName)
                    if (length(which(data[,w+2]==0)) > 0) {
                        data[,w+2][which(data[,w+2]==0)] = 0.001
                        wellsWithZeros = c(wellsWithZeros, wellName)
                        message(sprintf("***  %-3s contains zero-value data points  ***", wellName))
                    }
                    if (length(which(is.na(data[,w+2]))) > 0) {
                        data[,w+2][which(is.na(data[,w+2]))] = 0.0001
                        wellsWithZeros = c(wellsWithZeros, wellName)
                        message(sprintf("***  %-3s contains no data                 ***", wellName))
                    }
                    smooth = rbind(smooth, data.frame(Time=data$Time, OD=data[,w+2], SmoothedOD=avgSmooth(data[,w+2], smoothWindowSize,log=F), Well=wellName))
                }
            }

            # calculate growth rates (1st derivative)
            message("\nRates (1st derivative)...")
#            rate = tapply(smooth$SmoothedOD, smooth$Well, ratesInterval, smooth$Time, 11, TRUE)
            rate = tapply(smooth$SmoothedOD, smooth$Well, ratesInterval, smooth$Time, 15, TRUE)
            r = cbind(smooth, Rate=unlist(rate), NRate=unlist(rate))

            # calculate 2nd derivative
            message("\nRates' (2nd derivative)...")
#            rate1 = tapply(r$Rate, r$Well, ratesInterval, data$Time, 9, FALSE)
            rate1 = tapply(r$Rate, r$Well, ratesInterval, data$Time, 15, FALSE)
            r1 = cbind(r, NRate1=unlist(rate1))


            # calculate 3rd derivative
            message("\nRates'' (3rd derivative)...")
#            rate2 = tapply(r1$NRate1, r1$Well, ratesOffsets, data$Time, 3, 5, FALSE)
            rate2 = tapply(r1$NRate1, r1$Well, ratesOffsets, data$Time, 4, 6, FALSE)
            newdata <<- cbind(r1, NRate2=unlist(rate2))


            wellNames <<- as.vector(unique(newdata$Well))
        } else {

            ############### transformations specific for "08202009-2.txt" data set from Selwyn ###############

            tmp=read.table(dataFileName, header=T)
            data = data.frame(Time=sapply(tmp$Time, seconds), OD=tmp$OD, Well=tmp$Well)
            data$Well=as.vector(data$Well)
            times=data$Time[data$Well=="A01"]
            data = rbind(data, createMissingWell("A05", times, 0.01))
            data = rbind(data, createMissingWell("A10", times, 0.01))
            data = rbind(data, createMissingWell("A11", times, 0.01))
            data = rbind(data, createMissingWell("A12", times, 0.01))
            data = rbind(data, createMissingWell("B05", times, 0.01))
            data = rbind(data, createMissingWell("B10", times, 0.01))
            data = rbind(data, createMissingWell("B11", times, 0.01))
            data = rbind(data, createMissingWell("B12", times, 0.01))
            data = rbind(data, createMissingWell("C05", times, 0.01))
            data = rbind(data, createMissingWell("C10", times, 0.01))
            data = rbind(data, createMissingWell("C11", times, 0.01))
            data = rbind(data, createMissingWell("C12", times, 0.01))
            data = rbind(data, createMissingWell("D05", times, 0.01))
            data = rbind(data, createMissingWell("D10", times, 0.01))
            data = rbind(data, createMissingWell("D11", times, 0.01))
            data = rbind(data, createMissingWell("D12", times, 0.01))
            data = rbind(data, createMissingWell("E05", times, 0.01))
            data = rbind(data, createMissingWell("E10", times, 0.01))
            data = rbind(data, createMissingWell("E11", times, 0.01))
            data = rbind(data, createMissingWell("E12", times, 0.01))
            data = rbind(data, createMissingWell("F05", times, 0.01))
            data = rbind(data, createMissingWell("F10", times, 0.01))
            data = rbind(data, createMissingWell("F11", times, 0.01))
            data = rbind(data, createMissingWell("F12", times, 0.01))
            data = rbind(data, createMissingWell("G01", times, 0.01))
            data = rbind(data, createMissingWell("G02", times, 0.01))
            data = rbind(data, createMissingWell("G03", times, 0.01))
            data = rbind(data, createMissingWell("G04", times, 0.01))
            data = rbind(data, createMissingWell("G05", times, 0.01))
            data = rbind(data, createMissingWell("G06", times, 0.01))
            data = rbind(data, createMissingWell("G07", times, 0.01))
            data = rbind(data, createMissingWell("G08", times, 0.01))
            data = rbind(data, createMissingWell("G09", times, 0.01))
            data = rbind(data, createMissingWell("G10", times, 0.01))
            data = rbind(data, createMissingWell("G11", times, 0.01))
            data = rbind(data, createMissingWell("G12", times, 0.01))
            data = rbind(data, createMissingWell("H01", times, 0.01))
            data = rbind(data, createMissingWell("H02", times, 0.01))
            data = rbind(data, createMissingWell("H03", times, 0.01))
            data = rbind(data, createMissingWell("H04", times, 0.01))
            data = rbind(data, createMissingWell("H05", times, 0.01))
            data = rbind(data, createMissingWell("H06", times, 0.01))
            data = rbind(data, createMissingWell("H07", times, 0.01))
            data = rbind(data, createMissingWell("H08", times, 0.01))
            data = rbind(data, createMissingWell("H09", times, 0.01))
            data = rbind(data, createMissingWell("H10", times, 0.01))
            data = rbind(data, createMissingWell("H11", times, 0.01))
            data = rbind(data, createMissingWell("H12", times, 0.01))

            data = data[order(as.vector(sapply(data$Well, threeCharWell)), data$Time),]

            Rate=tapply(data$OD, data$Well, ratesInterval, times, 9)
            r=cbind(data, Rate=unlist(Rate))

            message("\nRates'...")
            rate1=tapply(r$Rate, r$Well, ratesInterval, times, 7, FALSE)
            r1=cbind(r, NRate1=unlist(rate1))

            message("\nRates''...")
            rate2=tapply(r1$NRate1, r1$Well, ratesOffsets, times, 1, 2, FALSE)
            newdata = cbind(r1, NRate2=unlist(rate2))

            wellNames <<- unique(newdata$Well)
        }


        message("\nAnalyzing wells...")
        analysisFileHeader = c("Well", "LagTime1", "GrowthRate1", "LagTime2", "LagTime2Ext", "GrowthRate2")
        analysisFileName = paste(dataFile,"_analysis",".txt",sep="")
        if (is.null(singleWell)) {
            write(analysisFileHeader, analysisFileName, ncolumns=6)
        }

        for (wellName in wellNames) {
            if (wellName %in% wellsWithZeros) {
                message(wellName, ":\t***  contains zero-value data points: not analyzed      ***")
                if (is.null(singleWell)) {
                    write(c(wellName, "NA", "NA", "NA", "NA"), analysisFileName, ncolumns=5, append=T)
                }
            } else {
                newdata$NRate[newdata$Well==wellName] <<- normalize(newdata$NRate[newdata$Well==wellName])
                newdata$NRate1[newdata$Well==wellName] <<- normalize(newdata$NRate1[newdata$Well==wellName])
                newdata$NRate2[newdata$Well==wellName] <<- normalize(newdata$NRate2[newdata$Well==wellName])

                # check if growth occured in well
                nTimePoints = length(newdata$OD[newdata$Well==wellName])
                meanODLastQuarter = mean(newdata$OD[newdata$Well==wellName][floor(nTimePoints*3/4):nTimePoints])
                meanODFirstTenth = mean(newdata$OD[newdata$Well==wellName][3:floor(nTimePoints*1/10)])
                hasGrowth = meanODLastQuarter >= (meanODFirstTenth * 1.333)

                if (!hasGrowth) {
                    message(wellName, ":\t***  no or less than 33% growth detected: not analyzed  ***")
                    if (is.null(singleWell)) {
                        write(c(wellName, "NA", "NA", "NA", "NA"), analysisFileName, ncolumns=5, append=T)
                    }
                } else {
                    message(wellName, ":")
                    wellProblem = FALSE
                    if (abs(log(mean(newdata$OD[newdata$Well==wellName][1:2]) / mean(newdata$OD[newdata$Well==wellName][3:5]))) > log(1.05)) {
                        wellProblem = TRUE
                        baseODdataPoints = 3:5
                    } else {
                        baseODdataPoints = 2:5
                    }

                    baseOD = mean(log(newdata$OD[newdata$Well==wellName][baseODdataPoints]))
                    baseODs <<- rbind(baseODs, data.frame(Well=wellName, OD=baseOD))

                    twoResourcesDetected = FALSE
                    if (isTwoResourceWell(oneResourceWells, twoResourceWells, wellName)) {
                        # analysis for wells with 2 resources
                        tmp = newdata[newdata$Well==wellName,]

                        valleyRanks = rankValley(tmp$Time, tmp$SmoothedOD, tmp$NRate, tmp$NRate1, tmp$NRate2)
                        tmpV = cbind(tmp, Rank=valleyRanks)
                        tmpV = tmpV[order(-tmpV$Rank),]
                        valley = tmpV[1,]
                        valleyIndex = which(tmp$Time==valley$Time[1])

                        tmpLeft = tmp[tmp$Time<=valley$Time[1],]
                        peak1Ranks = rankPeaks(tmpLeft$OD, tmpLeft$NRate, tmpLeft$NRate1, tmpLeft$NRate2)
                        tmpP1 = cbind(tmpLeft, Rank=peak1Ranks)
                        # exclude first 7 data points
                        tmpP1 = tmpP1[8:length(tmpP1$Time),]
                        tmpP1 = tmpP1[order(-tmpP1$Rank),]
                        peak1 = tmpP1[1,]
                        peak1Index = which(tmp$Time==peak1$Time[1])

                        expGrowth1 = which(tmpLeft$Rate > peak1$Rate[1]-peak1$Rate[1]/6)
                        expGrowth1 = expGrowth1[expGrowth1 > 2]

                        # cut data points from exponential growth interval
                        # if left or right side of interval has considerably
                        # more data points than other side
                        expGrowth1Left = length(which(expGrowth1 < peak1Index))
                        expGrowth1Right = length(which(expGrowth1 > peak1Index))
                        # message("\t\t", expGrowth1Left, " : ", expGrowth1Right)
                        if (expGrowth1Left >= (expGrowth1Right * 2)) {
                            # message("\t\t", sprintf("%i ", expGrowth1))
                            expGrowth1 = expGrowth1[(1 + (expGrowth1Left - expGrowth1Right) %/% 2):length(expGrowth1)]
                            # message("\t\t-->")
                            # message("\t\t", sprintf("%i ", expGrowth1))
                        } else {
                            if (expGrowth1Right >= (expGrowth1Left * 2)) {
                                # message("\t\t", sprintf("%i ", expGrowth1))
                                expGrowth1 = expGrowth1[1:(length(expGrowth1) - (expGrowth1Right - expGrowth1Left) %/% 2)]
                                # message("\t\t-->")
                                # message("\t\t", sprintf("%i ", expGrowth1))
                            }
                        }
                        peak1Time = minutesAndSeconds(peak1$Time[1])
                        message(sprintf("\t%-18s %5imin %02.0fs  (index: %3i, # data points: %2i)", "growth phase 1 at:", peak1Time[1], peak1Time[2], peak1Index, length(expGrowth1)))
                        message("\t\t", sprintf("%i ", expGrowth1))

                        if (peak1$Rate[1]*60*60 < 0.5) {
                            peak1MinWidth = 9
                        } else {
                            peak1MinWidth = 5
                        }
                        if (is.na(expGrowth1)) {
                            message("\t\t-->")
                            expGrowth1 = (which(tmpLeft$Rate==peak1$Rate[1])-peak1MinWidth%/%2) : (which(tmpLeft$Rate==peak1$Rate[1])+peak1MinWidth%/%2)
                            message("\t\t", sprintf("%i ", expGrowth1))
                        }
                        # add data points to both sides of exponential growth interval
                        # if length of interval is smaller than peak1MinWidth data points
                        expGrowth1Length = length(expGrowth1)
                        if (expGrowth1Length < peak1MinWidth) {
                            message("\t\t-->")
                            expGrowthAdd = (peak1MinWidth - expGrowth1Length + 1) %/% 2
                            expGrowth1 = c(seq(expGrowth1[1]-expGrowthAdd, expGrowth1[1]-1), expGrowth1, seq(expGrowth1[expGrowth1Length]+1, expGrowth1[expGrowth1Length]+expGrowthAdd))
                            message("\t\t", sprintf("%i ", expGrowth1))
                        }
                        expGrowth1 = expGrowth1[expGrowth1 > 2]


                        valleyTime = minutesAndSeconds(valley$Time[1])
                        message(sprintf("\t%-18s %5imin %02.0fs  (index: %3i)", "valley at:", valleyTime[1], valleyTime[2], valleyIndex))


                        tmpRight = tmp[tmp$Time>=valley$Time[1],]
                        peak2Ranks = rankPeaks(tmpRight$OD, tmpRight$NRate, tmpRight$NRate1, tmpRight$NRate2)
                        tmpP2 = cbind(tmpRight, Rank=peak2Ranks)
                        tmpP2 = tmpP2[order(-tmpP2$Rank),]
                        peak2 = tmpP2[1,]
                        peak2Index = which(tmp$Time==peak2$Time[1])

                        expGrowth2 = which(tmpRight$Rate > peak2$Rate[1]-peak2$Rate[1]/6) + which(tmp$Time==valley$Time[1]) - 1

                        # cut data points from exponential growth interval
                        # if left or right side of interval has considerably
                        # more data points than other side
                        expGrowth2Left = length(which(expGrowth2 < peak2Index))
                        expGrowth2Right = length(which(expGrowth2 > peak2Index))
                        # message("\t\t", expGrowth2Left, " : ", expGrowth2Right)
                        if (expGrowth2Left >= (expGrowth2Right * 2)) {
                            # message("\t\t", sprintf("%i ", expGrowth2))
                            expGrowth2 = expGrowth2[(1 + (expGrowth2Left - expGrowth2Right) %/% 2):length(expGrowth2)]
                            # message("\t\t-->")
                            # message("\t\t", sprintf("%i ", expGrowth2))
                        } else {
                            if (expGrowth2Right >= (expGrowth2Left * 2)) {
                                # message("\t\t", sprintf("%i ", expGrowth2))
                                expGrowth2 = expGrowth2[1:(length(expGrowth2) - (expGrowth2Right - expGrowth2Left) %/% 2)]
                                # message("\t\t-->")
                                # message("\t\t", sprintf("%i ", expGrowth2))
                            }
                        }
                        peak2Time = minutesAndSeconds(peak2$Time[1])
                        message(sprintf("\t%-18s %5imin %02.0fs  (index: %3i, # data points: %2i)", "growth phase 2 at:", peak2Time[1], peak2Time[2], peak2Index, length(expGrowth2)))
                        message("\t\t", sprintf("%i ", expGrowth2))

                        if (peak2$Rate[1]*60*60 < 0.5) {
                            peak2MinWidth = 9
                        } else {
                            if (peak2$Rate[1]*60*60 < 0.95) {
                                peak2MinWidth = 5
                            } else {
                                peak2MinWidth = 3
                            }
                        }
                        if (is.na(expGrowth2)) {
                            message("\t\t-->")
                            expGrowth2 = ((which(tmpRight$Rate==peak2$Rate[1])-peak2MinWidth%/%2) : (which(tmpRight$Rate==peak2$Rate[1])+peak2MinWidth%/%2)) + which(tmp$Time==valley$Time[1]) - 1
                            message("\t\t", sprintf("%i ", expGrowth2))
                        }
                        # add data points to both sides of exponential growth interval
                        # if length of interval is smaller than peak1MinWidth data points
                        expGrowth2Length = length(expGrowth2)
                        if (expGrowth2Length < peak2MinWidth) {
                            message("\t\t-->")
                            expGrowthAdd = (peak2MinWidth - expGrowth2Length + 1) %/% 2
                            expGrowth2 = c(seq(expGrowth2[1]-expGrowthAdd, expGrowth2[1]-1), expGrowth2, seq(expGrowth2[expGrowth2Length]+1, expGrowth2[expGrowth2Length]+expGrowthAdd))
                            message("\t\t", sprintf("%i ", expGrowth2))
                        }


                        # test if 2 separate growth phases were detected
                        if ((abs(valleyIndex - peak1Index) <= 5) || (abs(valleyIndex - peak2Index) <= 5) || ((valley$OD[1] > peak2$OD[1]) && (valley$NRate[1] < -0.15))) {
                            message("\n\t*****  UNABLE TO DETECT SEPARATE GROWTH PHASES FOR 2 RESOURCES  *****")
                            if ((abs(valleyIndex - peak1Index) <= 5) || (abs(valleyIndex - peak2Index) <= 5)) {
                                message("\t*****  valley and peak too close together                       *****\n")
                            } else {
                                message("\t*****  OD of valley > OD of peak 2 and high negative rate       *****\n")
                            }
                        } else {
                            twoResourcesDetected = TRUE

                            ranksValley <<- rbind(ranksValley, data.frame(Well=wellName, Time=tmp$Time, Rank=valleyRanks))
                            ranksPeak1 <<- rbind(ranksPeak1, data.frame(Well=wellName, Time=tmpLeft$Time, Rank=peak1Ranks))
                            ranksPeak2 <<- rbind(ranksPeak2, data.frame(Well=wellName, Time=tmpRight$Time, Rank=peak2Ranks))

                            valleys <<- rbind(valleys, valley)
                            peak1s <<- rbind(peak1s, peak1)
                            peak2s <<- rbind(peak2s, peak2)
                            expGrowths <<- rbind(expGrowths, data.frame(Well=wellName, Left1=expGrowth1[1], Peak1=peak1Index, Right1=expGrowth1[length(expGrowth1)], Left2=expGrowth2[1], Peak2=peak2Index, Right2=expGrowth2[length(expGrowth2)]))

                            model1 = lm(log(newdata$OD[newdata$Well==wellName][expGrowth1]) ~ newdata$Time[newdata$Well==wellName][expGrowth1])
                            # print(summary(model))
                            r21 = summary(model1)$r.squared
#                             message("\t", r21)
                            c1 = coef(model1)

                            model2 = lm(log(newdata$OD[newdata$Well==wellName][expGrowth2]) ~ newdata$Time[newdata$Well==wellName][expGrowth2])
                            # print(summary(model))
                            r22 = summary(model2)$r.squared
#                             message("\t", r22)
                            c2 = coef(model2)

                            models[[length(models)+1]] <<- list(Model1=model1, Model2=model2)
                            names(models) <<- c(names(models)[1:(length(models)-1)], wellName)


                            lag1StartTime = newdata$Time[newdata$Well==wellName][1]
                            lag1EndTime = (baseOD-c1[1])/c1[2]
                            lag1 = lag1EndTime - lag1StartTime

                            lag2OD = mean(log(tmp$OD[(valleyIndex-2):(valleyIndex+2)]))
                            lag2EndTime = (lag2OD-c2[1])/c2[2]
                            lag2 = lag2EndTime - valley$Time[1]

                            lag2StartIndex = valleyIndex
                            for (i in (valleyIndex - 1):peak1Index) {
                                if ((tmp$SmoothedOD[i] - tmp$SmoothedOD[lag2StartIndex]) >= -0.0003) {
                                    lag2StartIndex = i;
                                } else {
                                    break;
                                }
                            }
                            lag2StartTimeExtended = tmp$Time[lag2StartIndex]
                            lag2extended = lag2EndTime - lag2StartTimeExtended

                            lags <<- rbind(lags, data.frame(Well=wellName, BaseOD=baseOD, Lag1StartTime=lag1StartTime, Lag1EndTime=lag1EndTime, Lag2OD=lag2OD, Lag2StartTime=valley$Time[1], Lag2StartTimeExt=lag2StartTimeExtended, Lag2EndTime=lag2EndTime))

                            lag1mns = minutesAndSeconds(lag1)
                            lag2mns = minutesAndSeconds(lag2)
                            lag2extendedmns = minutesAndSeconds(lag2extended)
                            message("\n\t**")
                            message(sprintf("\t%-13s %9.3fs  (%imin %02.0fs)", "LagTime1:", lag1, lag1mns[1], lag1mns[2]))
                            message(sprintf("\t%-13s %11.9f per hour", "GrowthRate1:", c1[2]*60*60))
                            message(sprintf("\t%-13s %9.3fs  (%imin %02.0fs)", "LagTime2:", lag2, lag2mns[1], lag2mns[2]))
                            message(sprintf("\t%-13s %9.3fs  (%imin %02.0fs)", "LagTime2Ext:", lag2extended, lag2extendedmns[1], lag2extendedmns[2]))
                            message(sprintf("\t%-13s %11.9f per hour", "GrowthRate2:", c2[2]*60*60))
                            if (is.null(singleWell)) {
                                write(c(wellName, sprintf("%.3f", lag1), sprintf("%.9f", c1[2]*60*60), sprintf("%.3f", lag2), sprintf("%.3f", lag2extended), sprintf("%.9f", c2[2]*60*60)), analysisFileName, ncolumns=6, append=T)
                            }

                            if (wellProblem) {
                                message("\n\t***  High initial drop or jump in OD value  ***")
                            }
                            if (r21 < rSquaredTreshold) {
                                if (!wellProblem) {
                                    message()
                                }
                                rSquared = formatC(r21, digits=3, flag="#", format = "f")
                                message("\t***  Weak fit of regression line for growth phase 1 (R^2 = ", rSquared, " < ", rSquaredTreshold, ")  ***")
                                wellProblem = TRUE
                            }
                            if (r22 < rSquaredTreshold) {
                                if (!wellProblem) {
                                    message()
                                }
                                rSquared = formatC(r22, digits=3, flag="#", format = "f")
                                message("\t***  Weak fit of regression line for growth phase 2 (R^2 = ", rSquared, " < ", rSquaredTreshold, ")  ***")
                                wellProblem = TRUE
                            }
                            if (wellProblem) {
                                problems <<- rbind(problems, data.frame(Well=wellName))
                            }
                        }
                        message()
                    }


                    if (!twoResourcesDetected) {
                        # analysis for wells with 1 resource
                        tmp = newdata[newdata$Well==wellName,]
                        # exclude first 7 data points
                        tmp2 = tmp[8:length(tmp$Time),]
                        # maximum rate that has no strong negative rate'
                        maxRate = max(tmp2$Rate[tmp2$NRate1 > -0.2])
                        peakIndex = which(tmp$Rate == maxRate)
                        expGrowth = which(tmp$Rate > maxRate-maxRate/6)
                        # expGrowth = which(rates>quantile(rates[rates>mean(rates)],.7))
                        expGrowth = expGrowth[expGrowth > 2]

                        if (length(expGrowth) > 7) {
                            # remove data points from exponential growth interval
                            # that are too far off to the left or right of peak rate
                            expGrowthLeft = expGrowth[expGrowth < peakIndex]
                            expGrowthRight = expGrowth[expGrowth > peakIndex]
                            expGrowthLeftLength = length(expGrowthLeft)
                            expGrowthRightLength = length(expGrowthRight)
                            distance = c(peakIndex - expGrowthLeft, 0, expGrowthRight - peakIndex)
                            # print(distance)
                            expectedDistance = c(seq(expGrowthLeftLength,1), 0, seq(1,expGrowthRightLength))
                            # print(expectedDistance)
                            expGrowth = expGrowth[distance - expectedDistance < 5]
                        }

                        # cut data points from exponential growth interval
                        # if left or right side of interval has considerably
                        # more data points than other side
                        expGrowthLeftLength = length(expGrowth[expGrowth < peakIndex])
                        expGrowthRightLength = length(expGrowth[expGrowth > peakIndex])
                        # message("\t\t", expGrowthLeftLength, " : ", expGrowthRightLength)
                        if (expGrowthLeftLength >= (expGrowthRightLength * 1.8)) {
                            # message("\t\t", sprintf("%i ", expGrowth))
                            expGrowth = expGrowth[(1 + (expGrowthLeftLength - expGrowthRightLength) %/% 1.2):length(expGrowth)]
                            # message("\t\t-->")
                            # message("\t\t", sprintf("%i ", expGrowth))
                        } else {
                            if (expGrowthRightLength >= (expGrowthLeftLength * 1.8)) {
                                # message("\t\t", sprintf("%i ", expGrowth))
                                expGrowth = expGrowth[1:(length(expGrowth) - (expGrowthRightLength - expGrowthLeftLength) %/% 1.2)]
                                # message("\t\t-->")
                                # message("\t\t", sprintf("%i ", expGrowth))
                            }
                        }

                        # add data points to both sides of exponential growth interval
                        # if length of interval is smaller than 9 data points
                        expGrowthLength = length(expGrowth)
                        if (expGrowthLength < 9) {
                            expGrowthAdd = (9 - expGrowthLength + 1) %/% 2
                            expGrowth = c(seq(expGrowth[1]-expGrowthAdd, expGrowth[1]-1), expGrowth, seq(expGrowth[expGrowthLength]+1, expGrowth[expGrowthLength]+expGrowthAdd))
                        }
                        expGrowth = expGrowth[expGrowth > 2]

                        model = lm(log(tmp$OD[expGrowth]) ~ tmp$Time[expGrowth])
#                         print(summary(model))
                        r2 = summary(model)$r.squared
#                         message("\t", r2)
                        c = coef(model)

                        expGrowths <<- rbind(expGrowths, data.frame(Well=wellName, Left1=expGrowth[1], Peak1=peakIndex, Right1=expGrowth[length(expGrowth)], Left2=NA, Peak2=NA, Right2=NA))
                        models[[length(models)+1]] <<- list(Model1=model)
                        names(models) <<- c(names(models)[1:(length(models)-1)], wellName)

                        lagStartTime = tmp$Time[1]
                        lagEndTime = (baseOD-c[1])/c[2]
                        lag = lagEndTime - lagStartTime
                        lags <<- rbind(lags, data.frame(Well=wellName, BaseOD=baseOD, Lag1StartTime=lagStartTime, Lag1EndTime=lagEndTime, Lag2OD=NA, Lag2StartTime=NA, Lag2StartTimeExt=NA, Lag2EndTime=NA))

                        lagmns = minutesAndSeconds(lag)
                        message(sprintf("\t%-13s %9.3fs  (%imin %02.0fs)", "LagTime:", lag, lagmns[1], lagmns[2]))
                        message(sprintf("\t%-13s %11.9f per hour", "GrowthRate:", c[2]*60*60))
                        if (is.null(singleWell)) {
                            write(c(wellName, sprintf("%.3f", lag), sprintf("%.9f", c[2]*60*60), "NA", "NA", "NA"), analysisFileName, ncolumns=6, append=T)
                        }

                        if (wellProblem) {
                            message("\n\t***  High initial drop or jump in OD value  ***")
                        }
                        if (r2 < rSquaredTreshold) {
                            if (!wellProblem) {
                                message()
                            }
                            rSquared = formatC(r2, digits=3, flag="#", format = "f")
                            message("\t***  Weak fit of regression line (R^2 = ", rSquared, " < ", rSquaredTreshold, ")  ***")
                            wellProblem = TRUE
                        }
                        if (wellProblem) {
                            problems <<- rbind(problems, data.frame(Well=wellName))
                        }
                        message()
                    }
                }
            }
        }
    }



    message("\nGenerating figure(s):")

    if (is.null(singleWell)) {
        pdfHeight = 8
        pdfWidth = 15
        pdfWellName = ""
        grid = c(12,8)
        cexMultiplier = 1
        cexPointMultiplier = 1
        lwdMultiplier = 1
    } else {
        pdfHeight = 5
        pdfWidth = 6
        pdfWellName = paste("_single_", singleWell, sep="")
        grid = c(1,1)
        cexMultiplier = 4.5
        cexPointMultiplier = cexMultiplier/3
        lwdMultiplier = 1.5
    }

    minRate = min(newdata$Rate[!newdata$Well %in% wellsWithZeros])*60*60
    maxRate = max(newdata$Rate[!newdata$Well %in% wellsWithZeros])*60*60
    ratesYLim = c(minRate-(maxRate-minRate)*0.075, maxRate+(maxRate-minRate)*0.075)
    ratesYTicks = c(-0.4, 0, 0.4, 0.8)

    if (plotDetail) {
        pdfName = paste(dataFile, pdfWellName, "_detail.pdf", sep="")
    } else {
        pdfName = paste(dataFile, pdfWellName, ".pdf", sep="")
    }
    pdf(pdfName, height=pdfHeight, width=pdfWidth)
    par(mar=c(4,5,2,1)+0.2, cex.lab=1.4, cex.axis=1.2, cex.main=1)


    message("\nGenerating ODs figure...")

    panel.index <- 1

    plot = xyplot(OD ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab="Time (h)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(log="e", at=yTicks, tck=0.7)), panel = function(x, y, subscripts,...) {

        wellName = wellNames[panel.index]

        # check if growth occured in well
        hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
        hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

        if (!hasGrowth) {
            if (wellName %in% wellsWithZeros) {
                # plot red background for wells with problems
                panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,0,0,0.12))
            } else {
                # plot yellow background for wells with no growth
                panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,1,0,0.12))
            }
        } else {
            baseOD = baseODs$OD[baseODs$Well==wellName]

            # plot dashed green line for base OD level
            panel.segments(minX, baseOD, maxX*0.5, baseOD, lty=2, col=modelCol, lwd=lwdMultiplier*0.666)

            if (hasTwoGrowthPhases) {
                # plot dashed green line for OD level of beginning of 2nd lag time (valley)
                panel.segments(lags$Lag2StartTimeExt[lags$Well==wellName]-maxX*0.05, lags$Lag2OD[lags$Well==wellName], lags$Lag2EndTime[lags$Well==wellName]+maxX*0.15, lags$Lag2OD[lags$Well==wellName], lty=2, col=modelCol, lwd=lwdMultiplier*0.666)

                # plot red background for wells with problems
                if (length(problems$Well[problems$Well==wellName]) != 0) {
                    panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,0,0,0.12))
                }

                if (plotDetail) {
                    panel.abline(v=x[expGrowths$Left1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.3)
                    panel.abline(v=x[expGrowths$Right1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.3)
                    panel.abline(v=x[expGrowths$Left2[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.3)
                    panel.abline(v=x[expGrowths$Right2[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.3)

                    panel.points(valleys$Time[valleys$Well==wellName], log(valleys$OD[valleys$Well==wellName]), col="red", cex=cexPointMultiplier*0.6, lwd=lwdMultiplier*0.5)
                    panel.points(peak1s$Time[peak1s$Well==wellName], log(peak1s$OD[peak1s$Well==wellName]), col=modelCol, cex=cexPointMultiplier*0.6, lwd=lwdMultiplier*0.5)
                    panel.points(peak2s$Time[peak2s$Well==wellName], log(peak2s$OD[peak2s$Well==wellName]), col=modelCol, cex=cexPointMultiplier*0.6, lwd=lwdMultiplier*0.5)
                }
            } else {
                # plot blue background for 2-resource wells with single growth phase
                if (isTwoResourceWell(oneResourceWells, twoResourceWells, wellName)) {
                    panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(0,0,1,0.06))
                }

                # plot red background for wells with problems
                if (length(problems$Well[problems$Well==wellName]) != 0) {
                    panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,0,0,0.12))
                }

                if (plotDetail) {
                    panel.abline(v=x[expGrowths$Left1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.3)
                    panel.abline(v=x[expGrowths$Right1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.3)

                    panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, cex=cexPointMultiplier*0.6, lwd=lwdMultiplier*0.5)
                }
            }
        }

        # plot well name
        panel.text(x=minX, y=log(wellNameY), pos=4, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
        # plot 2 large dots under well name for wells with 2 resources
        if (isTwoResourceWell(oneResourceWells, twoResourceWells, wellName)) {
            panel.points(minX+(maxX-minX)*0.075, log(twoResourceY), pch=19, cex=cexMultiplier*0.5, col="gray40", lwd=lwdMultiplier)
            if (hasTwoGrowthPhases) {
                panel.points(minX+(maxX-minX)*0.135, log(twoResourceY), pch=19, cex=cexMultiplier*0.5, col="gray40", lwd=lwdMultiplier)
            } else {
                panel.points(minX+(maxX-minX)*0.135, log(twoResourceY), pch=21, cex=cexMultiplier*0.5, col="gray40", bg="white", lwd=lwdMultiplier*1.2)
            }
        }

        # plot growth data
        panel.xyplot(x, y, lwd=lwdMultiplier, ...)

        # plot models and estimates
        if (hasGrowth) {
            baseOD = baseODs$OD[baseODs$Well==wellName]

            # plot 1st lag time line
            lag1Y = baseOD + (quantile(y,0.9) - baseOD) * 0.38
            panel.arrows(lags$Lag1StartTime[lags$Well==wellName], lag1Y, lags$Lag1EndTime[lags$Well==wellName], lag1Y, angle=90, length=lwdMultiplier*0.02, code=3, lwd=lwdMultiplier*0.8, lend=2, col=modelCol)

            if (hasTwoGrowthPhases) {
                # plot 2nd lag time line
                lag2Y = lags$Lag2OD[lags$Well==wellName] + (quantile(y,0.95) - lags$Lag2OD[lags$Well==wellName]) * 0.98
                panel.arrows(lags$Lag2StartTimeExt[lags$Well==wellName], lag2Y, lags$Lag2StartTime[lags$Well==wellName], lag2Y, angle=90, length=lwdMultiplier*0.02, code=3, lwd=lwdMultiplier*0.666, lend=2, col=modelColLight)
                panel.arrows(lags$Lag2StartTime[lags$Well==wellName], lag2Y, lags$Lag2EndTime[lags$Well==wellName], lag2Y, angle=90, length=lwdMultiplier*0.02, code=3, lwd=lwdMultiplier*0.8, lend=2, col=modelCol)

                c1 = coef(models[[wellName]][["Model1"]])
                c2 = coef(models[[wellName]][["Model2"]])

                # plot growth phase models lines
                panel.abline(c1, col=modelCol, lwd=lwdMultiplier*0.666)
                panel.abline(c2, col=modelCol, lwd=lwdMultiplier*0.666)

                # plot result estimates
                lag1mns = minutesAndSeconds(lags$Lag1EndTime[lags$Well==wellName] - lags$Lag1StartTime[lags$Well==wellName])
                lag2mns = minutesAndSeconds(lags$Lag2EndTime[lags$Well==wellName] - lags$Lag2StartTime[lags$Well==wellName])
                panel.text(x=maxX+5, y=log(lagValue1Y), pos=2, labels=sprintf("%im%02.0fs", lag1mns[1], lag1mns[2]), col="black", cex=cexMultiplier*0.47)
                panel.text(x=maxX+5, y=log(lagValue2Y), pos=2, labels=sprintf("%im%02.0fs", lag2mns[1], lag2mns[2]), col="black", cex=cexMultiplier*0.47)
                panel.text(x=maxX+5, y=log(slopeValue1Y), pos=2, labels=formatC(c1[2]*60*60, digits=5, flag="#", format = "f"), col="black", cex=cexMultiplier*0.47)
                panel.text(x=maxX+5, y=log(slopeValue2Y), pos=2, labels=formatC(c2[2]*60*60, digits=5, flag="#", format = "f"), col="black", cex=cexMultiplier*0.47)
            } else {

                c = coef(models[[wellName]][["Model1"]])

                # plot growth phase model line
                panel.abline(c, col=modelCol, lwd=lwdMultiplier*0.666)

                # plot result estimates
                lagmns = minutesAndSeconds(lags$Lag1EndTime[lags$Well==wellName] - lags$Lag1StartTime[lags$Well==wellName])
                panel.text(x=maxX+5, y=log(lagValueY), pos=2, labels=sprintf("%im%02.0fs", lagmns[1], lagmns[2]), col="black", cex=cexMultiplier*0.47)
                panel.text(x=maxX+5, y=log(slopeValueY), pos=2, labels=formatC(c[2]*60*60, digits=5, flag="#", format = "f"), col="black", cex=cexMultiplier*0.47)
            }
        }

        panel.index <<- panel.index + 1
    })

    print(plot)


    # plot derivatives
    if (plotDetail) {

        message("\nGenerating smoothed ODs figure...")

        panel.index <- 1

        plot = xyplot(SmoothedOD ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab="Time (h)", ylab="OD (smoothed)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(log="e", at=yTicks, tck=0.7)), panel = function(x, y, subscripts,...) {

            wellName = wellNames[panel.index]

            # check if growth occured in well
            hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
            hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

            if (!hasGrowth) {
                if (wellName %in% wellsWithZeros) {
                    # plot red background for wells with problems
                    panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,0,0,0.12))
                } else {
                    # plot yellow background for wells with no growth
                    panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,1,0,0.12))
                }
            } else {
                if (hasTwoGrowthPhases) {
                    # plot red background for wells with problems
                    if (length(problems$Well[problems$Well==wellName]) != 0) {
                        panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,0,0.0,0.12))
                    }

                    panel.points(valleys$Time[valleys$Well==wellName], log(valleys$SmoothedOD[valleys$Well==wellName]), col="red", cex=cexPointMultiplier*0.6)
                    panel.points(peak1s$Time[peak1s$Well==wellName], log(peak1s$SmoothedOD[peak1s$Well==wellName]), col=modelCol, cex=cexPointMultiplier*0.6)
                    panel.points(peak2s$Time[peak2s$Well==wellName], log(peak2s$SmoothedOD[peak2s$Well==wellName]), col=modelCol, cex=cexPointMultiplier*0.6)
                } else {
                    # plot blue background for 2-resource wells with single growth phase
                    if (isTwoResourceWell(oneResourceWells, twoResourceWells, wellName)) {
                        panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(0,0,1,0.06))
                    }

                    # plot red background for wells with problems
                    if (length(problems$Well[problems$Well==wellName]) != 0) {
                        panel.rect(minX, log(minY), maxX, log(maxY), col=rgb(1,0,0.0,0.12))
                    }

                    panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, cex=cexPointMultiplier*0.6)
                }
            }

            # plot well name
            panel.text(x=minX, y=log(wellNameY), pos=4, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
            # plot 2 large dots under well name for wells with 2 resources
            if (isTwoResourceWell(oneResourceWells, twoResourceWells, wellName)) {
                panel.points(minX+(maxX-minX)*0.075, log(twoResourceY), pch=19, cex=cexMultiplier*0.5, col="gray40", lwd=lwdMultiplier)
                if (hasTwoGrowthPhases) {
                    panel.points(minX+(maxX-minX)*0.135, log(twoResourceY), pch=19, cex=cexMultiplier*0.5, col="gray40", lwd=lwdMultiplier)
                } else {
                    panel.points(minX+(maxX-minX)*0.135, log(twoResourceY), pch=21, cex=cexMultiplier*0.5, col="gray40", bg="white", lwd=lwdMultiplier*1.2)
                }
            }

            # plot growth data
            panel.xyplot(x, y, lwd=lwdMultiplier, ...)

            # plot models and estimates
            if (hasGrowth) {

                # plot 1st lag time line
                panel.arrows(lags$Lag1StartTime[lags$Well==wellName], lags$BaseOD[lags$Well==wellName], lags$Lag1EndTime[lags$Well==wellName], lags$BaseOD[lags$Well==wellName], angle=90, length=lwdMultiplier*0.02, code=3, lwd=lwdMultiplier*0.8, lend=2, col=modelCol)

                if (hasTwoGrowthPhases) {
                    # plot 2nd lag time line
                    panel.arrows(lags$Lag2StartTimeExt[lags$Well==wellName], lags$Lag2OD[lags$Well==wellName], lags$Lag2StartTime[lags$Well==wellName], lags$Lag2OD[lags$Well==wellName], angle=90, length=lwdMultiplier*0.02, code=3, lwd=lwdMultiplier*0.666, lend=2, col=modelColLight)
                    panel.arrows(lags$Lag2StartTime[lags$Well==wellName], lags$Lag2OD[lags$Well==wellName], lags$Lag2EndTime[lags$Well==wellName], lags$Lag2OD[lags$Well==wellName], angle=90, length=lwdMultiplier*0.02, code=3, lwd=lwdMultiplier*0.8, lend=2, col=modelCol)

                    panel.abline(v=x[expGrowths$Left1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)
                    panel.abline(v=x[expGrowths$Right1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)
                    panel.abline(v=x[expGrowths$Left2[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)
                    panel.abline(v=x[expGrowths$Right2[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)

                    # plot growth phase models lines
                    panel.abline(coef(models[[wellName]][["Model1"]]), col=modelCol, lwd=lwdMultiplier*0.666)
                    panel.abline(coef(models[[wellName]][["Model2"]]), col=modelCol, lwd=lwdMultiplier*0.666)

                    panel.text(x=maxX+5, y=log(lagValue1Y), pos=2, labels=paste("p:", formatC(coef(summary(models[[wellName]][["Model1"]]))[2,4], digits=1, flag="#", format = "e")), col="black", cex=cexMultiplier*0.47)
                    panel.text(x=maxX+5, y=log(lagValue2Y), pos=2, labels=paste("p:", formatC(coef(summary(models[[wellName]][["Model2"]]))[2,4], digits=1, flag="#", format = "e")), col="black", cex=cexMultiplier*0.47)
                    rSquared = formatC(summary(models[[wellName]][["Model1"]])$r.squared, digits=3, flag="#", format = "f")
                    if (rSquared < rSquaredTreshold) {
                        rSquaredCol = "red2"
                    } else {
                        rSquaredCol = "black"
                    }
                    panel.text(x=maxX+5, y=log(slopeValue1Y), pos=2, labels=bquote(R^2: .(rSquared)), col=rSquaredCol, cex=cexMultiplier*0.47)
                    rSquared = formatC(summary(models[[wellName]][["Model2"]])$r.squared, digits=3, flag="#", format = "f")
                    if (rSquared < rSquaredTreshold) {
                        rSquaredCol = "red2"
                    } else {
                        rSquaredCol = "black"
                    }
                    panel.text(x=maxX+5, y=log(slopeValue2Y), pos=2, labels=bquote(R^2: .(rSquared)), col=rSquaredCol, cex=cexMultiplier*0.47)

                } else {
                    panel.abline(v=x[expGrowths$Left1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)
                    panel.abline(v=x[expGrowths$Right1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)

                    # plot growth phase model line
                    panel.abline(coef(models[[wellName]][["Model1"]]), col=modelCol, lwd=lwdMultiplier*0.666)

                    panel.text(x=maxX+5, y=log(slopeValue1Y), pos=2, labels=paste("p:", formatC(coef(summary(models[[wellName]][["Model1"]]))[2,4], digits=1, flag="#", format = "e")), col="black", cex=cexMultiplier*0.47)
                    rSquared = formatC(summary(models[[wellName]][["Model1"]])$r.squared, digits=3, flag="#", format = "f")
                    if (rSquared < rSquaredTreshold) {
                        rSquaredCol = "red2"
                    } else {
                        rSquaredCol = "black"
                    }
                    panel.text(x=maxX+5, y=log(slopeValue2Y), pos=2, labels=bquote(R^2: .(rSquared)), col=rSquaredCol, cex=cexMultiplier*0.47)
                }
            }

            panel.index <<- panel.index + 1
        })

        print(plot, newpage=TRUE)


        message("\nGenerating rates figure (1st derivative)...")

        panel.index <- 1

        plot = xyplot(Rate*60*60 ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=ratesYLim, xlab="Time (h)", ylab="Growth rate (per hour)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(tick.number=3, tck=0.7)), panel = function(x, y, subscripts, ...) {
        # plot = xyplot(Rate*60*60 ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=ratesYLim, xlab="Time (h)", ylab="Growth rate (per hour)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(at=ratesYTicks, tck=0.7)), panel = function(x, y, subscripts, ...) {

            wellName = wellNames[panel.index]
            hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

            panel.text(x=maxX+5, y=ratesYLim[1]+0.9*(ratesYLim[2]-ratesYLim[1]), pos=2, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
            panel.abline(h=0, lty=5, col="gray", lwd=lwdMultiplier*0.8)
            panel.xyplot(x, y, lwd=lwdMultiplier, ...)

            hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
            if (hasGrowth) {
                if (hasTwoGrowthPhases) {
                    panel.points(valleys$Time[valleys$Well==wellName], valleys$Rate[valleys$Well==wellName]*60*60, col="red", pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(valleys$Time[valleys$Well==wellName], valleys$Rate[valleys$Well==wellName]*60*60, col="red", cex=cexPointMultiplier*0.7)

                    panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$Rate[peak1s$Well==wellName]*60*60, col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$Rate[peak1s$Well==wellName]*60*60, col=modelCol, cex=cexPointMultiplier*0.7)
                    panel.abline(h=peak1s$Rate[peak1s$Well==wellName]*60*60-peak1s$Rate[peak1s$Well==wellName]*60*60/6, lty=3, col="red", lwd=lwdMultiplier*0.8)

                    panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$Rate[peak2s$Well==wellName]*60*60, col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$Rate[peak2s$Well==wellName]*60*60, col=modelCol, cex=cexPointMultiplier*0.7)
                    panel.abline(h=peak2s$Rate[peak2s$Well==wellName]*60*60-peak2s$Rate[peak2s$Well==wellName]*60*60/6, lty=3, col="red", lwd=lwdMultiplier*0.8)

                } else {
                    panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, cex=cexPointMultiplier*0.7)

                    # panel.abline(h=max(y)-0.001, lty=3, col="green", lwd=lwdMultiplier*0.8)
                    panel.abline(h=y[expGrowths$Peak1[expGrowths$Well==wellName]]-y[expGrowths$Peak1[expGrowths$Well==wellName]]/6, lty=3, col=modelCol, lwd=lwdMultiplier*0.8)
                    # panel.abline(h=quantile(y[y>mean(y)],.7), lty=3, col="red", lwd=lwdMultiplier*0.8)

                    panel.abline(v=x[expGrowths$Left1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)
                    panel.abline(v=x[expGrowths$Right1[expGrowths$Well==wellName]], lty=5, col=modelColLight, lwd=lwdMultiplier*0.4)
                }
            }

            panel.index <<- panel.index + 1
        })

        print(plot, newpage=TRUE)


        if (twoResourceWells != "NONE") {
            message("\nGenerating normalized rates figure (1st derivative)...")

            panel.index <- 1

            plot = xyplot(NRate ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=c(-0.7,1.1), xlab="Time (h)", ylab="Growth rate (normalized)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(at=c(-0.4, 0, 0.4, 0.8), tck=0.7)), panel = function(x, y, subscripts, ...) {

                wellName = wellNames[panel.index]
                hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

                panel.text(x=maxX+5, y=-0.7+0.9*(1.1+0.7), pos=2, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
                panel.abline(h=0, lty=5, col="gray", lwd=lwdMultiplier*0.8)
                panel.xyplot(x, y, lwd=lwdMultiplier, ...)

                hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
                if (hasGrowth) {
                    if (hasTwoGrowthPhases) {
                        panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate[valleys$Well==wellName], col="red", pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate[valleys$Well==wellName], col="red", cex=cexPointMultiplier*0.7)
                        panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate[peak1s$Well==wellName], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate[peak1s$Well==wellName], col=modelCol, cex=cexPointMultiplier*0.7)
                        panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate[peak2s$Well==wellName], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate[peak2s$Well==wellName], col=modelCol, cex=cexPointMultiplier*0.7)
                    } else {
                        panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, cex=cexPointMultiplier*0.7)
                    }
                }

                panel.index <<- panel.index + 1
            })

            print(plot, newpage=TRUE)
        }


        message("\nGenerating normalized rates' figure (2nd derivative)...")

        panel.index <- 1

        plot = xyplot(NRate1 ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=c(-0.9,0.9), xlab="Time (h)", ylab="Growth Rate' (normalized)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(at=c(-0.6, 0, 0.6), tck=0.7)), panel = function(x, y, subscripts, ...) {

            wellName = wellNames[panel.index]
            hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

            panel.text(x=maxX+5, y=-0.9+0.9*(0.9+0.9), pos=2, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
            panel.abline(h=0, lty=5, col="gray", lwd=lwdMultiplier*0.8)
            panel.xyplot(x, y, lwd=lwdMultiplier, ...)

            hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
            if (hasGrowth) {
                if (hasTwoGrowthPhases) {
                    panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate1[valleys$Well==wellName], col="red", pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate1[valleys$Well==wellName], col="red", cex=cexPointMultiplier*0.7)
                    panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate1[peak1s$Well==wellName], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate1[peak1s$Well==wellName], col=modelCol, cex=cexPointMultiplier*0.7)
                    panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate1[peak2s$Well==wellName], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate1[peak2s$Well==wellName], col=modelCol, cex=cexPointMultiplier*0.7)
                } else {
                    panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                    panel.points(x[expGrowths$Peak1[expGrowths$Well==wellName]], y[expGrowths$Peak1[expGrowths$Well==wellName]], col=modelCol, cex=cexPointMultiplier*0.7)
                }
            }

            panel.index <<- panel.index + 1
        })

        print(plot, newpage=TRUE)



        if (toupper(twoResourceWells) != "NONE") {
            message("\nGenerating normalized rates'' figure (3rd derivative)...")

            panel.index <- 1

            plot = xyplot(NRate2 ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=c(-0.9,0.9), xlab="Time (h)", ylab="Growth rate'' (normalized)", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(at=c(-0.6, 0, 0.6), tck=0.7)), panel = function(x, y, subscripts, ...) {

                wellName = wellNames[panel.index]
                hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

                panel.text(x=maxX+5, y=-0.9+0.9*(0.9+0.9), pos=2, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
                panel.abline(h=0, lty=5, col="gray", lwd=lwdMultiplier*0.8)
                panel.xyplot(x, y, lwd=lwdMultiplier, ...)

                hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
                if (hasGrowth) {
                    if (hasTwoGrowthPhases) {
                        panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate2[valleys$Well==wellName], col="red", pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate2[valleys$Well==wellName], col="red", cex=cexPointMultiplier*0.7)
                        panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate2[peak1s$Well==wellName], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate2[peak1s$Well==wellName], col=modelCol, cex=cexPointMultiplier*0.7)
                        panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate2[peak2s$Well==wellName], col=modelCol, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate2[peak2s$Well==wellName], col=modelCol, cex=cexPointMultiplier*0.7)
                    }
                }

                panel.index <<- panel.index + 1
            })

            print(plot, newpage=TRUE)



            message("\nGenerating ranks figure...")

            panel.index <- 1

            plot = xyplot(NRate ~ Time|Well, newdata, layout=grid, type='l', as.table=T, strip=F, xlim=c(minX,maxX), ylim=c(-2.2,1.2), xlab="Time (h)", ylab="Ranks of valley and peaks", main=paste(dataFile,".txt",sep=""), scales = list(x=list(at=xTicks, labels=timeTicks, tck=0.7), y=list(tick.number=0, tck=0.7)), panel = function(x, y, subscripts, ...) {

                wellName = wellNames[panel.index]
                hasTwoGrowthPhases = isTwoResourceWell(oneResourceWells, twoResourceWells, wellName) && (length(valleys$Well[valleys$Well==wellName]) != 0)

                panel.text(x=maxX+5, y=-2.2+0.9*(1.2+2.2), pos=2, labels=wellName, col="darkgray", cex=cexMultiplier*0.7)
                panel.abline(h=0, lty=5, col="lightgray", lwd=lwdMultiplier*0.6)
                panel.xyplot(x, y, lwd=lwdMultiplier*0.75, col="lightblue", ...)

                hasGrowth = (length(baseODs$Well[baseODs$Well==wellName]) != 0)
                if (hasGrowth) {
                    if (hasTwoGrowthPhases) {
                        panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate[valleys$Well==wellName], col="lightpink", pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(valleys$Time[valleys$Well==wellName], valleys$NRate[valleys$Well==wellName], col="lightpink", cex=cexPointMultiplier*0.7)
                        panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate[peak1s$Well==wellName], col=modelColLight, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(peak1s$Time[peak1s$Well==wellName], peak1s$NRate[peak1s$Well==wellName], col=modelColLight, cex=cexPointMultiplier*0.7)
                        panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate[peak2s$Well==wellName], col=modelColLight, pch=19, cex=cexPointMultiplier*0.08)
                        panel.points(peak2s$Time[peak2s$Well==wellName], peak2s$NRate[peak2s$Well==wellName], col=modelColLight, cex=cexPointMultiplier*0.7)

                        valleyRanks = ranksValley$Rank[ranksValley$Well==wellName]
                        valleyRanks = (valleyRanks - min(valleyRanks)) / (max(valleyRanks)-min(valleyRanks))
                        panel.lines(ranksValley$Time[ranksValley$Well==wellName], valleyRanks-2, lwd=lwdMultiplier, col="red")
                        peak1Ranks = ranksPeak1$Rank[ranksPeak1$Well==wellName]
                        peak1Ranks = (peak1Ranks - min(peak1Ranks)) / (max(peak1Ranks)-min(peak1Ranks))
                        panel.lines(ranksPeak1$Time[ranksPeak1$Well==wellName], peak1Ranks-1.4, lwd=lwdMultiplier, col=modelCol)
                        peak2Ranks = ranksPeak2$Rank[ranksPeak2$Well==wellName]
                        peak2Ranks = (peak2Ranks - min(peak2Ranks)) / (max(peak2Ranks)-min(peak2Ranks))
                        panel.lines(ranksPeak2$Time[ranksPeak2$Well==wellName], peak2Ranks-1.4, lwd=lwdMultiplier, col=modelCol)
                    }
                }

                panel.index <<- panel.index + 1
            })

            print(plot, newpage=TRUE)
        }
    }

    dev.off()

    runTime = proc.time() - ptm
    message("\nDone\nRun time: ", sprintf("%imin %02.0fs", runTime[3]%/%60, runTime[3]%%60))
}




# analyzeGrowthCurves("f2 and f8 clones in glu 28_7_02 linear.txt", twoResourceWells=c("B1","D1","B2","F2","B4","D4","D5","A6","B6","C6","D6","F6","G6","D7","B9","D9","F9","B10","D10","F10","C11","D11","F11","A12","B12","D12","E12","F12"), plotDetail=T)

# analyzeGrowthCurves("MEL-15_april.pda.txt", twoResourceWells=c(9, 11), plotDetail=T)


# analyzeGrowthCurves("06_19_11_Diauxic_Pause_14hrRun_REL606_LacY3.3-1.txt", twoResourceWells=c("ALL"), plotDetail=T, maxTimeHours=15, timeTicks=c(0,5,10,15))

# analyzeGrowthCurves("06_19_11_Diauxic_Pause_14hrRun_REL606_LacY3.3-1.txt", twoResourceWells=c("ALL"), plotDetail=T, maxTimeHours=15, timeTicks=c(0,5,10,15), maxWells=9)

# analyzeGrowthCurves("06_19_11_Diauxic_Pause_14hrRun_REL606_LacY3.3-1.txt", plotDetail=T)

# analyzeGrowthCurves("06_19_11_Diauxic_Pause_14hrRun_REL606_LacY3.3-1.txt")


# analyzeGrowthCurves("06_30_11_Diauxic_Pause_12hrRun_REL606.txt", twoResourceWells=c("ALL"), plotDetail=T, maxTimeHours=15, timeTicks=c(0,5,10,15))

# analyzeGrowthCurves("06_30_11_Diauxic_Pause_12hrRun_REL606.txt", twoResourceWells=c("ALL"), plotDetail=T)

# analyzeGrowthCurves("06_30_11_Diauxic_Pause_12hrRun_REL606.txt", twoResourceWells=c("ALL"), plotDetail=T, maxOD=0.3, ODTicks=c(0.1,0.2,0.3), maxTimeHours=12, timeTicks=c(0,4,8,12), singleWell="A7")

# analyzeGrowthCurves("06_30_11_Diauxic_Pause_12hrRun_REL606.txt", plotDetail=T, singleWell="A7")

# analyzeGrowthCurves("06_30_11_Diauxic_Pause_12hrRun_REL606.txt", plotDetail=T)

# analyzeGrowthCurves("06_30_11_Diauxic_Pause_12hrRun_REL606.txt")



analyzeGrowthCurves("VP_12h_DM_NaCl-Glu_Replicate2_10-31-11.txt", twoResourceWells=c("NONE"), maxOD=1.7, ODTicks=c(0.1,0.5,1.5), maxTimeHours=25, timeTicks=c(0,8,16,24), plotDetail=T)

