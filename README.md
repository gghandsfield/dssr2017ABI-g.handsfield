# Muscle Fiber Determination from Laplacian Simulations- A Comparison with DTI

## Conducting a Laplacian Simulation for Segmented Muscle Data
This step is performed using Autodesk Inventor CAD software  
Use Manage>Import to import STL muscle data into Autodesk Inventor  
Use the Tool to create a base feature from STL  
Create planes that transect the geometry at regular intervals  
In each plane, create 2D cubic splines as a 2D Sketch that approximate the geometry  
Note: you can use slice graphics to cut away part of the STL to see what you're doing  
When finished, loft the cubic splines to a solid

###### Creating Aponeuroses on your Muscle Model
You need to create regions of your muscle solid that correspond to the location of aponeuroses  
One approach to do this is to create a surface in 3D space that cuts off part of your solid  
The part of the solid that is cut off should be the aponeurosis part

###### Set up Model in Autodesk CFD
Use the icon on the 'Model' tab to send your model to Autodesk CFD  
Define proximal aponeurosis surface(s) as inlets with 1Pa gage pressure  
Define distal aponeurosis surface(s) as outlets with 0Pa gage pressure  
Define all other surface(s) as slip/symmetry  
Define material as 1g/cm^3, 1000Pa-s fluid  
Mesh the geometry using automatic feature  
Conduct simulation  
When complete, export nodal data in the File menu

## Loading Muscle Surfaces and Fibre Vectors from DTI
In order to perform comparisons, must have DTI eigenvector data and STL muscle surfaces from segmentation
