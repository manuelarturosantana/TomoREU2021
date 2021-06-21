# TomoREU2021: Point-of-care Tomographic Imaging
Code written to solve tomographic imaging problems during an REU/RET at Emory University during the summer of 2021.

The research team, led by Dr. James Nagy from Emory's Mathematics Department, consists of 4 student fellows:
1. Manuel Santana
2. Mai Phuong Pham Huynh
3. Issa Susa
4. Ana Castillo

In this project, we aim to use computational method to develop a numerical method to jointly estimate the geometry parameters of the portable CT scan device and to reconstruct the image.

For full functionality the following matlab packages and toolboxes should be installed.
* [Matlab Optmization Toolbox](https://www.mathworks.com/products/optimization.html)
* [IRtools](https://github.com/jnagy1/IRtools)
* [Imfill](https://ctk.math.ncsu.edu/imfil.html)

The GITHUB repository consists of 4 main folders:
1. test_image : image to set up test for our algorithm
2. test_setup : set up tests
3. algorithms : MATLAB functions
4. test_reslt : results produced from the test set-up (reconstructed images, relative error norm graphs)
