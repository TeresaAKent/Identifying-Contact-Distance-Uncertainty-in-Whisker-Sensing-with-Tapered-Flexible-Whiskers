# Identifying-Contact-Distance-Uncertainty-in-Whisker-Sensing-with-Tapered-Flexible-Whiskers
###Teresa A. Kent, Hannah Emnett, Mahnoush Babaei, Mitra J. Z. Hartmann, Sarah Bergbreiter

**Abstract-Whisker-based tactile sensors have the potential to perform fast and accurate 3D mappings of the environment, complementing vision-based methods under conditions of glare, reflection, proximity, and occlusion. Current algorithms for mapping with whiskers make assumptions about the conditions of contact, but these assumptions are not always valid and can cause significant sensing errors. Here we introduce a new whisker sensing system with a tapered, flexible whisker. The system provides inputs to two separate algorithms for estimating radial contact distance on a whisker. Using a Gradient-Moment (GM) algorithm, we correctly detect contact distance in most cases (within 4% of the whisker length). We introduce the Z-Dissimilarity score as a new metric that quantifies uncertainty in the radial contact distance estimate using both the GM algorithm and a Moment-Force (MF) algorithm that exploits the tapered whisker design. Combining the two algorithms ultimately results in contact distance estimates more robust than either algorithm alone.**

Files Available-

1. Code implementation of the tracking algorithm which converts a test video into a csv file of the motion of tracked points
2. Code to preform analysis on the components of the sensor including
    1. Output of predicted whisker bending from the quasistatic simulator developed at northwestern university in the Sense Lab
    2. Analysis on the sensitivity of acrylic cut springs to forces produced by whisker bending

*Kent, T. A., Emnett, H., Babaei, M., Hartmann, M. J., & Bergbreiter, S. (2023). Identifying Contact Distance Uncertainty in Whisker Sensing with Tapered, Flexible Whiskers. International Conference on Robotics and Automation*

Any questions about this code can be directed to tkent@andrew.cmu.edu.

If anyone would like a copy of the papers, video or data from this paper, the authors are happy to oblige a request by email. 

-------------------------------------------------------------------------
The purpose of each of the author written codes is as follows:

1. HeatMapThicknessvsWidthSpring.m In this matlab script the affect of design decisions for the spring suspension on the minimum force and minimum moment are calculated. These values are compared to the minimum sensitivity required to detect radial contact distance on a whisker 1 mm (black) and 10 mm (gray) apart. The final chosen design parameters for the springs is represented by a red star. The output of this code can be seen in the image below. 

![](https://github.com/TeresaAKent/Identifying-Contact-Distance-Uncertainty-in-Whisker-Sensing-with-Tapered-Flexible-Whiskers/blob/aea899dcb1025df002a8b39883981188dbb4d1a1/Design%20Decisions.png)

2. CustomSquareTrackingWithTemplateCorrection.py: This code takes a video from the sensor as input and outputs a csv file of the locations of the tracked points and videos which show the tracking. A pictoral representation of the tracking algorithm can be seen in the picture below. 
![](https://github.com/TeresaAKent/Identifying-Contact-Distance-Uncertainty-in-Whisker-Sensing-with-Tapered-Flexible-Whiskers/blob/aea899dcb1025df002a8b39883981188dbb4d1a1/VisualTracking.png)
