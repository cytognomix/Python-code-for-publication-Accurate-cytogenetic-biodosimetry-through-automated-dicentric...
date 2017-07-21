import os, shutil, math, itertools
import xml.etree.ElementTree as xmlTree
from array import array
import numpy as np

def demoIt(sampleFile):
    """This function is only for demo purpose"""
    images = readADCISampleXml2Images(sampleFile)
    filtersZScores(images)
    print sampleFile+" loaded, "+str(len(images))+" images."
    print ""
    print "Images sorted by Group Bin Distance, descending from good to poor:"
    scores = []
    sortIndex = []
    for image in images:
        scores.append(image.groupBinVectorDistance())
    sortIndex = sorted(range(len(scores)), key=lambda k:scores[k])
    for index in sortIndex:
        print images[index].getImageFileName()+", Group Bin Distance: "+str(scores[index])
    print ""
    print "Images sorted by Combined Z Score, weight [4,3,4,5,2,1], descending from good to poor:"
    scores = []
    sortIndex = []
    for image in images:
        scores.append(image.filterScore(array('d', [4.0,3.0,4.0,5.0,2.0,1.0])))
    sortIndex = sorted(range(len(scores)), key=lambda k:scores[k])
    for index in sortIndex:
        print images[index].getImageFileName()+", Combined Z Score: "+str(scores[index])
    print ""

class metaphase:
    """Data structure for metaphase image"""
	
    def __init__(self, name, startR, endR, index):
        self.name = name
        self.startR = startR
        self.endR = endR
        self.size = endR+1-startR
        self.index = index

        self.lengthByWidth = 0.0
        self.candidateDensity = 0.0
        self.finiteDiff = 0.0
        # Valid object: green + red + yellow contour
        self.countValid = 0
        # Segmented object: Valid object + blue contour
        self.countSegmented = 0
        # z scores
        self.zLengthWidth = 0.0
        self.zCandidateDensity = 0.0
        self.zFiniteDiff = 0.0
        self.zCountObj = 0.0
        self.zCountSegmented = 0.0
        self.zCountValid = 0.0
        # DC counts 
        self.dcs = array('i', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])
        # Chromosome areas, sorted, excluding noise (overlap class 12)
        self.areas = []
        # Result matrix list of list, 0-3 bounding box, 4-14 svm score
        self.results = []
    def showInfo(self):
        print self.name
        print self.startR
        print self.endR

    def getImageFileName(self):
        tokens = self.name.split('/')
        return tokens[-1]

    def filterScore(self, weights = array('d', [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])) :
        """Order: lengthWidth, candidateDensity, finiteDiff, obj, segmented,
classified. obj, segmented will be using abs value"""
        score = weights[0]*self.zLengthWidth
        score = score + weights[1]*self.zCandidateDensity
        score = score - weights[2]*self.zFiniteDiff
        score = score + weights[3]*math.fabs(self.zCountObj)
        score = score + weights[4]*math.fabs(self.zCountSegmented)
        score = score - weights[5]*self.zCountValid;
        return score

    def groupBinVector(self, female=True) :
        """3 Element vecotr, count of chromosomes in AB, C, D-G groups, by area"""
        # If we know the gender, we could use different bound for male/female
        boundAB_C = 0.029
        boundC_DG = 0.02
        groupBinVector = [0,0,0]
        totalArea = sum(self.areas)
        if totalArea == 0 :
            return groupBinVector;
        for area in self.areas :
            if float(area)/totalArea > boundAB_C :
                groupBinVector[0] = groupBinVector[0] + 1
            elif float(area)/totalArea > boundC_DG :
                groupBinVector[1] = groupBinVector[1] + 1
            else :
                groupBinVector[2] = groupBinVector[2] + 1
        return groupBinVector

    def groupBinVectorDistance(self, female=True) :
        """calculate the euclidean distance between area and standardArea"""
        if female :
            standardArea = [10, 16, 20]
        else :
            standardArea = [10, 15, 21]
        array1 = np.array(self.groupBinVector())
        array2 = np.array(standardArea)
        return np.linalg.norm(array1-array2)

    def areaDistributionDistance(self, female=True) :
        """calculate the euclidean distance bewteen area of standardArea directly
        instead of using group bin"""
        if female :
            standardArea = [0.041129832, 0.041129832, 0.039973571, 0.039973571,
                            0.032705649, 0.032705649, 0.031384209, 0.031384209,
                            0.029732408, 0.029732408, 0.028080608, 0.028080608,
                            0.026428807, 0.026428807, 0.025768087, 0.025768087,
                            0.023951107, 0.023951107, 0.022794846, 0.022794846,
                            0.022134126, 0.022134126, 0.022299306, 0.022299306,
                            0.021968946, 0.021968946, 0.018830525, 0.018830525,
                            0.017674265, 0.017674265, 0.016848365, 0.016848365,
                            0.014866204, 0.014866204, 0.013709944, 0.013709944,
                            0.013214404, 0.013214404, 0.010571523, 0.010571523,
                            0.009745623, 0.009745623, 0.008424182, 0.008424182,
                            0.007763462, 0.007763462]
        else :
            standardArea = [0.041813602, 0.041813602, 0.040638119, 0.040638119,
                            0.03324937, 0.03324937, 0.031905961, 0.031905961,
                            0.0302267, 0.0302267, 0.028547439, 0.028547439,
                            0.026868178, 0.026868178, 0.024349286, 0.024349286,
                            0.023173804, 0.023173804, 0.022502099, 0.022502099,
                            0.022670025, 0.022670025, 0.022334173, 0.022334173,
                            0.019143577, 0.019143577, 0.017968094, 0.017968094,
                            0.017128463, 0.017128463, 0.01511335, 0.01511335,
                            0.013937867, 0.013937867, 0.013434089, 0.013434089,
                            0.009907641, 0.009907641, 0.010747271, 0.010747271,
                            0.007892527, 0.007892527, 0.008564232, 0.008564232,
                            0.026196474, 0.009571788]
        standardArea = sorted(standardArea, None, None, True)

        areaDistribution = []
        totalArea = sum(self.areas)
        for area in self.areas :
            areaDistribution.append(float(area)/totalArea)

        while len(areaDistribution) < len(standardArea) :
            areaDistribution.append(0.0)
        while len(standardArea) < len(areaDistribution) :
            standardArea.append(0.0)

        array1 = np.array(standardArea)*100.0        # fraction to percentage 
        array2 = np.array(areaDistribution)*100.0
        return np.linalg.norm(array1-array2)

    def areaDistributionPValue(self, method='KS', female=True) :
        """ Available method: KS"""
        if female :
            standardArea = [0.041129832, 0.041129832, 0.039973571, 0.039973571,
                            0.032705649, 0.032705649, 0.031384209, 0.031384209,
                            0.029732408, 0.029732408, 0.028080608, 0.028080608,
                            0.026428807, 0.026428807, 0.025768087, 0.025768087,
                            0.023951107, 0.023951107, 0.022794846, 0.022794846,
                            0.022134126, 0.022134126, 0.022299306, 0.022299306,
                            0.021968946, 0.021968946, 0.018830525, 0.018830525,
                            0.017674265, 0.017674265, 0.016848365, 0.016848365,
                            0.014866204, 0.014866204, 0.013709944, 0.013709944,
                            0.013214404, 0.013214404, 0.010571523, 0.010571523,
                            0.009745623, 0.009745623, 0.008424182, 0.008424182,
                            0.007763462, 0.007763462]
        else :
            standardArea = [0.041813602, 0.041813602, 0.040638119, 0.040638119,
                            0.03324937, 0.03324937, 0.031905961, 0.031905961,
                            0.0302267, 0.0302267, 0.028547439, 0.028547439,
                            0.026868178, 0.026868178, 0.024349286, 0.024349286,
                            0.023173804, 0.023173804, 0.022502099, 0.022502099,
                            0.022670025, 0.022670025, 0.022334173, 0.022334173,
                            0.019143577, 0.019143577, 0.017968094, 0.017968094,
                            0.017128463, 0.017128463, 0.01511335, 0.01511335,
                            0.013937867, 0.013937867, 0.013434089, 0.013434089,
                            0.009907641, 0.009907641, 0.010747271, 0.010747271,
                            0.007892527, 0.007892527, 0.008564232, 0.008564232,
                            0.026196474, 0.009571788]
        areaDistribution = []
        totalArea = sum(self.areas)
        for area in self.areas :
            areaDistribution.append(float(area)/totalArea)
            
        if method.lower() == 'ks' :
            d, p, ne = KSTest.kstest(areaDistribution, standardArea)
            return p
        else :
            print 'Invalid method in areaDistributionPValue()'
            return 0

def readADCISampleXml2Images(fileName, fpFlag=0x7E):
    """Read xml-formatted ADCI sample file"""
    sampleTree = xmlTree.parse(fileName)
    rootEle = sampleTree.getroot()
    version = rootEle.attrib['version']
    # skip other fields, directly go to images
    imagesEle = rootEle.find('Images')
    imageCount = 0
    startLine = 0
    images = []
    for imageEle in imagesEle.findall('Image'):
        imageName = imageEle.find('FileName').text
        imageCount += 1
        resultMatrixEle = imageEle.find('ResultMatrix')
        rows = int(resultMatrixEle.get('rows'))
        image = metaphase(imageName, startLine, startLine+rows-1, imageCount)
        startLine = startLine+rows
        imageStats = imageEle.find('ImageStats').text.split(';')
        image.lengthByWidth = float(imageStats[0])
        image.candidateDensity = float(imageStats[1])
        image.finiteDiff = float(imageStats[2])
        image.countValid = float(imageStats[3])
        image.countSegmented = float(imageStats[4])
        for chromosomeEle in resultMatrixEle.findall('Row'):
            line = chromosomeEle.text.split(';')
            for i in range(len(line)):
                line[i] = int(line[i])
            result = line[1:5]
            if len(line) == 8:  # all svm results are the same
                result.extend(line[-1:]*11)
            else:               # there are 11 svm results in line
                result.extend(line[-11:])
            image.results.append(result)
            area = line[5]
            overlapClass = line[6]
            if int(overlapClass)!=12:
                image.areas.append(int(area))
            for i in range(11):
                dcCode = result[i+4]
                if dcCode>0:
                    if (dcCode&fpFlag)==0:
                        if (dcCode&0x00001)==1:
                            image.dcs[i] = image.dcs[i]+1
        image.areas.sort(cmp=None, key=None, reverse=True)
        images.append(image)
    return images

def filtersZScores(images):
    """calculate the z score of all 6 filters, for a given set of images
    (all images or filtered images in a sample)"""
    data = np.empty([len(images), 6], dtype=float)
    for index, image in enumerate(images) :
        data[index, 0] = image.lengthByWidth
        data[index, 1] = image.candidateDensity
        data[index, 2] = image.finiteDiff
        data[index, 3] = float(image.size)
        data[index, 4] = float(image.countSegmented)
        data[index, 5] = float(image.countValid)
    means = np.mean(data, axis=0, dtype=np.float64)
    sds = np.std(data, axis=0, dtype=np.float64)
    for index, image in enumerate(images) :
        image.zLengthWidth = (data[index, 0]-means[0])/sds[0]
        image.zCandidateDensity = (data[index, 1]-means[1])/sds[1]
        image.zFiniteDiff = (data[index, 2]-means[2])/sds[2]
        image.zCountObj = (data[index, 3]-means[3])/sds[3]
        image.zCountSegmented = (data[index, 4]-means[4])/sds[4]
        image.zCountValid = (data[index, 5]-means[5])/sds[5]

if __name__=="__main__":
    demoIt("./DemoSample1.adcisample")
    demoIt("./DemoSample2.adcisample")
