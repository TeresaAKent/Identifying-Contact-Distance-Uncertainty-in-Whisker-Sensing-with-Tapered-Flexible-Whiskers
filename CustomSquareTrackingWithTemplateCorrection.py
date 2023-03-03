import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from LucasKanade import LucasKanade
from LucasKanade import checkRect
import cv2
import os
import skimage.morphology

def Updateframe(frame):
    # First white balance
    # Second color switch
    frame=frame-np.min(frame)
    WhiteBalance=255/np.average(frame[50:250,-250:-50,:],(0,1))*frame
    #print(np.average(frame[50:250,-250:-50,:],(0,1)))
    #frame2=(WhiteBalance[:,:,0]-np.average(WhiteBalance[:,:,0]))-(WhiteBalance[:,:,2]-np.average(WhiteBalance[:,:,2]))
    #frame2=((WhiteBalance[:,:,0]-np.average(WhiteBalance,2))-(WhiteBalance[:,:,2]-np.average(WhiteBalance,2)))/np.average(WhiteBalance,2)
    frame2=((WhiteBalance[:,:,0]-np.average(WhiteBalance[:,:,:],2))/(np.average(WhiteBalance[:,:,:],2)+.1)-(WhiteBalance[:,:,2]-np.average(WhiteBalance[:,:,:],2))/(np.average(WhiteBalance[:,:,:],2)+.2))
    frame2-=np.min(frame2)
    frame2=frame2*255/np.max(frame2)
    return frame2.astype('uint8')

def SelectBoxes(frame,bboxes):
    while True:

        frame3=frame.copy()
    
    
    
        LongText='\n Starting from the top left corner of the desired box, drag a box centered at one of the tracked points, \n when happy with the drawn box press the space button \n to select another box press the space bar a second time \n on the last box selection only press the space bar once then press q'
        ShortText='See Console: Pick One, Press Space'
        ShortText2='To Pick Another Press Space again, Done press Q'
        frame3=cv2.putText(frame3,ShortText,(50,50),cv2.FONT_HERSHEY_SIMPLEX, 1,(0,0,0),2,cv2.LINE_AA)
        frame3=cv2.putText(frame3,ShortText2,(50,100),cv2.FONT_HERSHEY_SIMPLEX, 1,(0,0,0),2,cv2.LINE_AA)
        print(LongText)
        #cv2.imshow('check',frame3)

        bbox = cv2.selectROI('MultiTracker', frame3,fromCenter=True) 
        bboxes.append(bbox)
        print (bboxes)
    
        k = cv2.waitKey(0) & 0xFF
        if (k == 113):  # q is pressed
            break
    print('Selected bounding boxes {}'.format(bboxes))
    
    return bboxes

def ArrowVisualization(Video1OutputName,fourcc,framesPerSecond,ArrowMultiplier,FinalData,FramesToVideo):
    numberOfBoxes=np.size(FinalData,1)//2

    
    video=cv2.VideoWriter(Video1OutputName,fourcc,framesPerSecond,(np.size(FramesToVideo,1),np.size(FramesToVideo,0)))
    InitalPositions=(FinalData[2,:]).flatten()
    for j in range(0,np.size(FramesToVideo,3)-1):
        ArrowFrame=FramesToVideo[:,:,:,j].copy()
        ArrowFrame=ArrowFrame-np.min(ArrowFrame)
        ArrowFrame=255/np.average(ArrowFrame[50:250,-250:-50,:],(0,1))*ArrowFrame
        ArrowFrame=np.clip(ArrowFrame,0,255)
        ArrowFrame=np.ascontiguousarray(ArrowFrame, dtype=np.uint8)
        if j==0:
            video.write(np.uint8(ArrowFrame))
        else:
            for k in range(numberOfBoxes):
                XCenter=int(FinalData[j,k*2])
                YCenter=int(FinalData[j,k*2+1])
                ArrowX=XCenter+ArrowMultiplier*int(XCenter-InitalPositions[k*2])
                ArrowY=YCenter+ArrowMultiplier*int(YCenter-InitalPositions[k*2+1])
                cv2.arrowedLine(ArrowFrame,(XCenter,YCenter),(ArrowX,ArrowY),(255,255,0),5,tipLength=.4)
            # cv2.imshow( "Display window", img);
            # cv2.waitKey(0)
            video.write(np.uint8(ArrowFrame))
            
        if j%50==0:
            print('Arrow Video Percent', j/np.size(FramesToVideo,3))
            
    
    video.release()
    cv2.destroyAllWindows()
    
def GetRodLocation(frame,margin,bboxes,params):
    # Preform whitebalencing
    frame=frame-np.min(frame,(0,1))
    FrameCC=255/np.average(frame[0:125,-250:-50,:],(0,1))*frame
    FrameCC=np.clip(FrameCC[225:350,700:850,:],0,255)
    #plt.figure('FrameCC')
    #plt.imshow(FrameCC.astype('uint8'))
    #plt.figure()
    #FrameCC=(WhiteMultiplier*frame[175:400,500:1000,:]).astype(int)
    #print(np.shape(FrameCC))
    BandWFrame=abs(FrameCC[:,:,0]-FrameCC[:,:,2])
    params["currentFrame"]=FrameCC
    # First figure out the groups
    selem1=skimage.morphology.disk(6)
    # Get the saturation and value
    WLower=(0,0,0)
    WUpper=(120,70,90)
    
    Step3 = cv2.inRange(FrameCC, WLower, WUpper)
    #cv2.imshow('MaskA.jpg', Step3)
    Step3BandW=cv2.inRange(BandWFrame, 0, 80)
    #cv2.imshow('MaskB.jpg', Step3BandW)
    Step3b=Step3*Step3BandW
    #cv2.imshow('MaskFull2.jpg', Step3b)
    Step4=skimage.morphology.binary_closing(Step3b,selem1)
    Step4=skimage.morphology.binary_closing(Step4,selem1)
    Step5,numLet=skimage.measure.label(Step4,return_num=True,background=0)
    #Sets the box size for the helper image
    bufRange=2*margin
    for region in skimage.measure.regionprops(Step5):
        if region.area>=25:
            # create a box to add to bboxes
            minr,minc,maxr,maxc=region.bbox
            bbox=(minc-margin+700,minr-margin+225,maxc-minc+bufRange,maxr-minr+bufRange)
            bboxes.append(bbox)
            # pick a color to label it whisker or dot

    
    # Get the index values for the whiskers
    
    if np.size(bboxes)==0:
        bboxes=params["last"]
    params["last"]=bboxes
    bboxes=np.asarray(bboxes)
    #print(bboxes)
    return bboxes  
 
def TrackerVerification(Video2OutputName,fourcc,framesPerSecond,FramesToVideo,FinalData):
    numberOfBoxes=np.size(FinalData,1)//2

        
    video=cv2.VideoWriter(Video2OutputName,fourcc,framesPerSecond,(np.size(FramesToVideo,1),np.size(FramesToVideo,0)))
    
    jUpperBound=int(np.size(FramesToVideo,3))
    for j in range(0,jUpperBound):
        img=FramesToVideo[:,:,:,j]
        img=np.ascontiguousarray(img, dtype=np.uint8)
        if j==0:
            video.write(np.uint8(img))
        else:
            for k in range(numberOfBoxes):
                XCenter=int(FinalData[j,k*2])
                YCenter=int(FinalData[j,k*2+1])
    
                #SizeX=FinalDataSizes[j,k*2]
                SizeX=40
                #SizeY=FinalDataSizes[j,k*2+1]
                SizeY=40
                
    
                
                p1 = (int(XCenter-SizeX/2), int(YCenter-SizeY/2))
                p2 = (int(XCenter+SizeX/2), int(YCenter+SizeY/2))
    
                cv2.rectangle(img,p1,p2,(0,0,0),2,1)
                #centroid_step.append(cent)
                img=cv2.circle(img,(XCenter,YCenter),10,(0,0,0))
                # print(int(i),int(cent[0]),int(cent[1]))
                TextString="{}".format(k)
                img=cv2.putText(img,TextString,(p2[0],p1[1]),cv2.FONT_HERSHEY_SIMPLEX, 1,(0,0,0),2,cv2.LINE_AA)
                
                video.write(np.uint8(img))
        if j%50==0:
            print('Tracker Verification Video Percent',j/np.size(FramesToVideo,3))
    video.release()
    cv2.destroyAllWindows()
    
#Load video
#Load first image  
videoPath = "../NSFTracker/Feburary 23rd/PhiTestTheta90MT11H32.5.mp4"
#Output folder name
OutputFolder="../NSFTracker/Feburary 23rd/BetterArrows/"
OutputfileName="PhiTestTheta90MT11H32.5"


cap = cv2.VideoCapture(videoPath)
framesPerSecond=cap.get(cv2.CAP_PROP_FPS)
TotalFrames=int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
print(TotalFrames)
success, frame = cap.read()
skipNumber=3
TotFrames=TotalFrames//skipNumber+2
FramesToVideo=np.zeros((np.size(frame,0),np.size(frame,1),3,TotFrames))

FramesToVideo[:,:,:,0]=frame
frame=Updateframe(frame)
frame0=frame.copy()
params={}


    
# Use select boxes to get a rectangle size
rectFull = np.asarray(SelectBoxes(frame,[])).T
oldrect=rectFull.copy()

maxx=np.size(frame,1)
maxy=np.size(frame,0)


carseqrects=np.zeros([TotFrames,4,np.size(oldrect,1)+1])
carseqrects[0,:,:-1]=rectFull
blackdot=GetRodLocation(FramesToVideo[:,:,:,0],25,[],params)
carseqrects[0,:,-1]=blackdot[-1,:]
imageNum=1
img0=frame
fig, ax = plt.subplots(nrows=1, ncols=1)

colors=['b','k','w','c','m','g','r','y']
hold=1
for i in range (TotalFrames-2):
    success, img=cap.read()
    if i%skipNumber==0:
        FramesToVideo[:,:,:,hold]=img
        blackdot=GetRodLocation(img,25,[],params)
        carseqrects[hold,:,-1]=blackdot[-1,:]
        img=Updateframe(img)
        
        for kk in range (np.size(oldrect,1)):
            rect=rectFull[:,kk]
            TempxSize=rect[2]
            TempySize=rect[3]
            if hold%TotalFrames ==500:
                # rounding
                rect[0]=int(round(rect[0]))
                rect[1]=int(round(rect[1]))
                imgBox=patches.Rectangle((rect[0],rect[1]),rect[2],rect[3], edgecolor=colors[kk], facecolor="none")
                imgBox2=patches.Rectangle((oldrect[0,kk],oldrect[1,kk]),oldrect[2,kk],oldrect[3,kk], edgecolor='cyan', facecolor="none")
                plt.imshow(img/np.max(img),cmap=plt.cm.gray)
                plt.gca().add_patch(imgBox)
                plt.gca().add_patch(imgBox2)
                imageNum=imageNum+1
                #add the current rectangle shape to it
            Temp=frame[rect[1]:rect[1]+rect[3],rect[0]:rect[0]+rect[2]]
            p0=np.zeros([2])
            PNew=LucasKanade(img0,img,rect)
            #PNew[0]=int(round(PNew[0]))
            #PNew[1]=int(round(PNew[1]))
            
            rect[0]=rect[0]+PNew[0]
            #rect[2]=int(round(rect[2]+PNew[0]))
            rect[1]=rect[1]+PNew[1]
            #rect[3]=int(round(rect[3]+PNew[1]))
            
            rect=checkRect(rect,img,TempxSize,TempySize)
            if i%5==0:
        
                PTempShift=LucasKanade(frame,img,oldrect[:,kk].flatten(),PNew)
            
                PNew[0]=PTempShift[0]-rect[0]+oldrect[0,kk]
                PNew[1]=PTempShift[1]-rect[1]+oldrect[1,kk]
                
                rect[0]=rect[0]+PNew[0]
                #rect[2]=int(round(rect[2]-PNew[0]))
                rect[1]=rect[1]+PNew[1]
                #rect[3]=int(round(rect[3]-PNew[1]))
                
                #rect=checkRect(rect,img,TempxSize,TempySize)
            
            rectFull[:,kk]=rect
            carseqrects[hold,:,kk]=rect
        print(i)
        img0=img
        hold+=1
    else:
        pass
#np.save("../code/girlseqrects-wcrt.npy",carseqrects)

# Arrow multiplier magnifies the movement of the dots for better visulization
ArrowMultiplier=3



if not os.path.exists(OutputFolder):
    os.mkdir(OutputFolder)
    print("Directory " , OutputFolder ,  " Created ")
else:    
    print("Directory " , OutputFolder,  " already exists")
OutputCSVFileName="{}{}Centers.csv".format(OutputFolder,OutputfileName)
Video1OutputName="{}{}ArrowVideo.mp4".format(OutputFolder,OutputfileName)
Video2OutputName="{}{}TrackingVerificationVideo.mp4".format(OutputFolder,OutputfileName)

# Output the csv file with point position

xvals=carseqrects[:,0,:]+carseqrects[:,2,:]/2
yvals=carseqrects[:,1,:]+carseqrects[:,3,:]/2
FinalData=np.zeros((np.size(xvals,0),np.size(xvals,1)+np.size(yvals,1)))
for i in range (np.size(FinalData,1)):
    if i%2==0:
        FinalData[:,i]=xvals[:,i//2]
    else:
        FinalData[:,i]=yvals[:,i//2]

np.savetxt(OutputCSVFileName, FinalData,delimiter=',')

fourcc = cv2.VideoWriter_fourcc(*'mp4v')

TrackerVerification(Video2OutputName,fourcc,framesPerSecond,FramesToVideo,FinalData)


ArrowVisualization(Video1OutputName,fourcc,framesPerSecond,ArrowMultiplier,FinalData,FramesToVideo)