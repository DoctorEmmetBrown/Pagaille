import numpy as np

from scipy import stats


def findCor(project0, project180, sizeBeamCor=20):
    shapeProjection = project180.shape
    shapeProjection = project0.shape

    maxPos = 0
    sizeBeam = sizeBeamCor; # Calculation of the COR on 10 lines

    intercorr = np.empty([sizeBeam, shapeProjection[1]])

    middle = int(shapeProjection[0] / 2)

    beg = int(middle - sizeBeamCor / 2)
    end = int(middle + sizeBeamCor / 2)
    for ii in range(beg, end): #np.arange(0,Itera,1):
        #print ii
        intercorr[ii - beg, :] = abs(np.correlate(project0[ii, :], project180[ii, :], mode='same', old_behavior=False))
        maxPos += intercorr[ii - beg].argmax()
        #print intercorr[ii-20].argmax()

    #plt.figure(0)
    #plt.plot(intercorr[5,:])
    #plt.show()
    #print intercorr[5].argmax()

    CORpos = maxPos / sizeBeam
    #    CORpos=1024
    print('COR found at  ' + str(CORpos))
    return CORpos


def findCorLudo(project0, project180, LineOfMaxEnergy, sizeBeamCor=20):
# shapeProjection=project180.shape
# shapeProjection=project0.shape

    maxPos = 0

    SumCorr = 0

    #middle=int(shapeProjection[0]/2)

    #print middle

    beg = int(LineOfMaxEnergy - sizeBeamCor / 2)
    end = int(LineOfMaxEnergy + sizeBeamCor / 2)

    Itera = end - beg
    intercorr = np.empty([Itera, project180.shape[1] * 2 - 1])
#    plt.figure(0)
#    plt.subplot(2, 1, 1)
    NbVox10pc = project180.shape[1] * 0.77

    corFoundsOnEachLine = np.arange(sizeBeamCor)

    for ii in np.arange(0, Itera, 1):

        prj180 = project180[ii + beg, :]
        prj0 = project0[ii + beg, :]
        hist180 = np.histogram(prj180, bins=400, range=None, normed=False, weights=None)
        cumSum180 = np.cumsum(hist180[0])
        valueThread180 = hist180[1][400 - len(cumSum180[cumSum180 > NbVox10pc])]
        #print 'valueThread180 '+str(valueThread180)
        #prj180[prj180 < valueThread180] = 0
        #prj180[prj180 >= valueThread180] = 1
        hist0 = np.histogram(prj0, bins=400, range=None, normed=False, weights=None)
        cumSum0 = np.cumsum(hist0[0])
        valueThread0 = hist0[1][400 - len(cumSum0[cumSum0 > NbVox10pc])]
        #prj0[prj0 < valueThread0] = 0
        #prj0[prj0 >= valueThread0] = 1
        #print 'valueThread0 '+str(valueThread0)

        #plt.plot(prj180)

        intercorr[ii, :] = np.correlate(prj0, prj180, mode='full')
        #plt.plot(intercorr[ii,:])
        SumCorr += np.max(intercorr[ii, :])
        #print(intercorr[ii].argmax() / 2.0)
        maxPos += (intercorr[ii].argmax()) * np.max(intercorr[ii, :])
        corFoundsOnEachLine[ii] = intercorr[ii].argmax() / 2.

    CORpos = np.true_divide(maxPos, SumCorr)
    CORpos /= 2.0
    CORpos = int(stats.mode(corFoundsOnEachLine)[0])
    #print "Find center of rotation at " + str(CORpos)
    return CORpos


def findCorNemoz(project0, project180, LineOfMaxEnergy, sizeBeamCor=20):
# shapeProjection=project180.shape
# shapeProjection=project0.shape

    maxPos = 0

    SumCorr = 0

    #middle=int(shapeProjection[0]/2)

    #print middle

    beg = int(LineOfMaxEnergy - sizeBeamCor / 2)
    end = int(LineOfMaxEnergy + sizeBeamCor / 2)

    Itera = end - beg
    intercorr = np.empty([Itera, project180.shape[1] * 2 - 1])
    plt.figure(0)
    NbVox10pc = project180.shape[1] * 0.9

    taMereLaMatrice = project0

    for ii in np.arange(0, Itera, 1):
        #print 'line '+str(ii)
        interCorSubMatrix = np.zeros((65, 4095))
        taMereLaMatrice = project0
        for gem in np.arange(16, 1024, 16):
            print('gem' + str(gem))
            prj180 = project180[ii + beg, :]
            prj0 = np.ones(2048)

            print('A un momnet va falloir qu on comprenne ' + str(prj0[prj0 != 0].size))

            prj0 = taMereLaMatrice[ii + beg, :]
            # test = prj0[prj0!=0]
            print('Au tout debut ' + str(taMereLaMatrice[taMereLaMatrice != 0].size))

            print('Alors que project0 fait ' + str(taMereLaMatrice.size))

            plt.plot(prj0)
            prj0[gem:2048] = 0

            print('Apres la feinte ' + str(prj0[prj0 != 0].shape))

            hist180 = np.histogram(prj180, bins=400, range=None, normed=False, weights=None)
            cumSum180 = np.cumsum(hist180[0])
            valueThread180 = hist180[1][400 - len(cumSum180[cumSum180 > NbVox10pc])]
            #print 'valueThread180 '+str(valueThread180)
            #prj180[prj180 < valueThread180] = 0
            #prj180[prj180 >= valueThread180] = 1
            hist0 = np.histogram(prj0, bins=400, range=None, normed=False, weights=None)
            cumSum0 = np.cumsum(hist0[0])
            valueThread0 = hist0[1][400 - len(cumSum0[cumSum0 > NbVox10pc])]
            #prj0[prj0 < valueThread0] = 0
            #prj0[prj0 >= valueThread0] = 1
            #print 'valueThread0 '+str(valueThread0)
            #plt.plot(prj0)
            #plt.plot(prj180)
            interCorSubMatrix[gem / 16, :] = np.correlate(prj0, prj180, mode='full', old_behavior=False)

            intercorr[ii, :] = np.correlate(prj0, prj180, mode='full', old_behavior=False)

            SumCorr += np.max(intercorr[ii, :])
            #  print intercorr[ii].argmax()/2.0
            maxPos += (intercorr[ii].argmax()) * np.max(intercorr[ii, :])
            print(intercorr[ii].argmax() / 2.)
        plt.show()
        # plt.imshow(interCorSubMatrix)
        #plt.show()

    CORpos = np.true_divide(maxPos, SumCorr)
    CORpos /= 2.0

    print( "Find center of rotation at " + str(CORpos))
    return CORpos

