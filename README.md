Random Field Topic Model
========

1, vc_RTM
It contains the Visual C++ source codes of the Random Topic Field Model used in our paper. 
The binary file is in release folder, along with the raw files of tracklets. The binary executable will first read the raw files of tracklets, and do hundreds of Gibbs Sampling iterations, which could be very slow. The output is classqq.txt. Please use the codes in the following matlab_IO folder to visualize it. 

2, matlab_IO
It contains the visualization tools for tracklets and semantic regions.

3, matlab_data
It contains the tracklet data in the format of matlab data. You could directly load it in the matlab. These tracklets are the same as the raw files of tracklets used in the release folder of vc_RTM. 

Note: 
1, the original video of the grandcentral station has been released at   		http://www.ee.cuhk.edu.hk/~xgwang/grandcentral.html

2, please cite our paper if you use any of the codes and data.
Bolei Zhou, Xiaogang Wang and Xiaoou Tang. "Random Field Topic Model for Semantic Region Analysis in Crowded Scenes from Tracklets." in Proceedings of IEEE Conference on Computer Vision and Pattern Recognition (CVPR 2011 )

Feel free to send me email at zhoubolei@gmail.com

Bolei 
Nov. 12, 2012
