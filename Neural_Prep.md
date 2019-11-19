Objective is to prepare a neural network to identify coordinates of a parallelogram which covers the dotted pattern in the image.

To define a parallelogram in 2D plane we need 6 parameters. Within our task the parameters selected to define the fig are.

1. Coordinates of opposite corners of paralellogram (X1, Y1) & (X3, Y3).
2. And length of sides of parallelogram L1 and L2.

As all images are taken at fixed distance and dotted pattern has fixed dimensions, we can reduce our parameters to 4 by keeping length of sides of parallelogram fixed.
Overall we need 4 coordinates (X1, Y1) & (X3, Y3) as output of neural network to define a parallelogram.

Dataset Prepration
To train the model we need to prepare dataset of approx. 500 images (in intital stage) of size 640 X 480 and a .csv file containing information of coordinates of opposite corners (X1, Y1) & (X3, Y3) of parrallelogram.

Model Prepration
For this task, we have selected 2 architectures
1. Vgg16
2. Resnet50

To train our model, we will use tranfer learning technique on both archetectures and select the model which is having better performance.

Objective is to prepare a neural network to identify coordinates of a parallelogram which covers the dotted pattern in the image.

To define a parallelogram in a 2D plane we need 6 parameters. Within our task the parameters selected to define the fig are.

1. Coordinates of opposite corners of a parallelogram (X1, Y1) & (X3, Y3).
2. And the length of sides of parallelogram L1 and L2.

As all images are taken at the fixed distance and dotted pattern has fixed dimensions, we can reduce our parameters to 4 by keeping the length of sides of the parallelogram fixed.
Overall we need 4 coordinates (X1, Y1) & (X3, Y3) as an output from the neural network to define a parallelogram.

Dataset Preparation
To train the model we need to prepare dataset of approx. 500 images (in initial stage) of size 640 X 480 and a .csv file containing information of coordinates of opposite corners (X1, Y1) & (X3, Y3) of parallelogram.

Model Preparation
For this task, we have selected 2 architectures
1. Vgg16
2. Resnet50
To train our model, we will use the transfer learning technique on both architectures and select the model which is having better performance.
