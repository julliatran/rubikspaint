from PIL import Image
import numpy as np
from scipy.misc import toimage, imread, imsave, imresize

def getFactors(x):
    factors = []
    for i in range(1, x+1):
        if x % i == 0: factors.append([i,x//i])
    return factors

#https://stackoverflow.com/questions/9018016/how-to-compare-two-colors-for-similarity-difference?fbclid=IwAR0iYGrU3vLo4sLFQf1MAUlM317CZcuarQYTMSAbBGKZma5BqLEQv-8Yji0
def getColorDistance(c0, c1):
    rmean = (c0[0] + c1[0]) / 2
    r = c0[0] - c1[0]
    g = c0[1] - c1[1]
    b = c0[2] - c1[2]
    return ((((512+rmean)*r*r)/256) + 4*g*g + (((767-rmean)*b*b)/256))**0.5

def imageToRubiks(filepath, maxNumCubes):
    imageData = imread(filepath)

    h, w = imageData.shape[0], imageData.shape[1]
    optimalRatio = float(w) / float(h)

    a, b = getFactors(int(maxNumCubes*0.8))[0]
    bestRatio = float(a) / float(b)
    bestRatioDifference = abs(bestRatio - optimalRatio)
    bestW = a
    bestH = b

    for numCubes in range(int(maxNumCubes*0.8)+1, maxNumCubes+1):
        factors = getFactors(numCubes)
        for a, b in factors:
            ratio = float(a) / float(b)
            ratioDifference = abs(ratio - optimalRatio)
            if ratioDifference < bestRatioDifference and a * b > bestW * bestH:
                bestRatio = ratio
                bestRatioDifference = ratioDifference
                bestW, bestH = a, b

    colorMatrix = imresize(imageData, (bestH * 3, bestW * 3), interp='bicubic')

    colors = [[196,30,58],[0,158,96],[0,81,186],
              [255,88,0],[255,213,0],[255,255,255]]

    for row in range(colorMatrix.shape[0]):
        for col in range(colorMatrix.shape[1]):
            pixelColor = colorMatrix[row,col]

            minDist = getColorDistance(pixelColor,colors[0])
            minIdx = 0

            for idx in range(1,len(colors)):
                colorDist = getColorDistance(pixelColor,colors[idx])
                if colorDist < minDist:
                    minDist = colorDist
                    minIdx = idx

            colorMatrix[row,col] = np.array(colors[minIdx])

    #previewData = imresize(colorMatrix, (12 * bestH, 12 * bestW), interp='nearest')
    previewImage = toimage(colorMatrix)
    previewFilepath = "%s_rubiks.png" % (filepath)
    previewImage.save(previewFilepath)

    return colorMatrix, previewFilepath
