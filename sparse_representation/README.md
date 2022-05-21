# 289FinalProject
The main references for this project are:

Yang, J., Wright, J., Huang, T. S., & Ma, Y. (2010). Image super-resolution via sparse representation. IEEE transactions on image processing, 19(11), 2861-2873.

Lee, H., Battle, A., Raina, R., & Ng, A. Y. (2007). Efficient sparse coding algorithms. In Advances in neural information processing systems (pp. 801-808).

## Learning Objective

While neural network sheds new insights in image processing in multiple applications including, image denoting, image blurring, image restoration, etc, conventional methods revoking linear algebra still remain active for 2 reasons, 1) linear algebra approach allows closer connection between theory and applications. 2) some machine learning ideas perfectly overlap and therefore empowers conventional image processing approaches. In this project, we are going to study one particular area of image reconstruction, namely image reconstruction to build a bridge between image processing and various types of least squares (least square regression problem, constrained least square problem), and demonstrate how to train a learning-based model through solving some linear programming problems. There are large amount of background knowledge in this note and we don’t expect the hypothetical students to understand them thoroughly. However, we aim to demonstrate a particular application of least squares method, specifically LASSO regression and Lagrangian Dual to enhance the students’ understanding of least square problem by an unfamiliar topic. We believe that applying what they learned in the class to unfamiliar problems is an important skill for the students to achieve a success career in industry. The background knowledge is provided to the students and the students will find out that the tasks left for them to do are what they have learned in previous weeks and 16B classes.

## Dictionaries
There are several important folders in the repo:
* ./Training: contains trained dictionary compact dictionary, and raw feature matrix.
* ./images: contains training images, demo example, and comparison.
* ./matlab_code: the LASSO regression and Constrained LS problems were solved with MATLAB scripts.
* ./python_code: the equivalent training code in python is given (it is not used in practice, it is provided for the sake of completeness and for the ease to demonstrate in notebook). The model is trained with MATLAB and the application of model is implemented use python.
* ./Report: contains final report of this project.
* ./slide: contains final PowerPoint and its pdf version.
* ./quiz: contains quiz or assignment.

## Content
## How to reproduce the results
The fully working jupyter notebook is in ./python_code/Sparse_Representation.ipynb. To run the code, downloads whole repo, otherwise the file address might not match. For some reasons, some Markdowns do not show in Github. We recommend users run it locally.


## Files in **matlab_code** folder:

* ./matlab_code/L1_FeatureSign.m：solver used to solve LASSO regression use feature sign algorithm.
* ./matlab_code/L1_FeatureSign_Setup.m: setup for L1_FeatureSign.m, call L1_FeatureSign.m for each column of sparse coding matrix separately, optimize it with L1_FeatureSign.m.
* ./matlab_code/L2_Lagrange_Dual.m: solver for Lagrange Dual problem for Constrained Least Square.
* ./matlab_code/objective_gradient_hessian.m: return objective function, gradient, and hessian.
* ./matlab_code/sparse_coding.m: set up for MATLAB optimization, input raw feature matrix. This function will return learned compact dictionary **D** and sparse coding **Z** by repeatedly calling L1_FeatureSign.m and L2_Lagrange_Dual.m to optimize **D** and **Z** iteratively.

## Files in **python_code** folder:
* ./python_code/Sparse_Representation.ipynb: jupyter-notebook contains everything, it is fully working. follow the procedure, one can get the same result showing in the notebook. The input image data and trained results are stored in folder images and Training, respectively.
* ./python_code/backprojection.py: used to further better performance by learning residual between input low resolution image and reconstructed high resolution image. This is optional.
* ./python_code/extract_low_resolution_features.py: extract first and second gradient in both horizontal and vertical direction from low resolution images.
* ./python_code/get_compact_X.py: read image from img_path, generate features and ensemble feature matrix.
* ./python_code/get_dictionary_demo.py: this is a demo code shows how to train and generate compact dictionary, it is not really useful, because the training process is completed in MATLAB.
* ./python_code/img_super_resolution.py: this is a wrap-up function used in application stage, it will return a high resolution image using learned compact dictionary from a low resolution image. this function will perform feature engineering for each patch (5x5 pixels) in low resolution image, run LASSO regression to find sparse coding for each patch, finally ensemble high resolution image.
* ./python_code/rnd_smp_patch.py: this function sets feature engineering steps for all images from the img_path, it will decide how many samples are taken from a image depending on its size and the total number of samples used to learn.
* ./python_code/sparse_coding.py: this function is the equivalent to the MATLAB scripts, it provides python version of two least squares problem solver. Unfortunately, it does not work properly because python does not have good linear programming solver for the proposed problem. If I have time, I will continue search for a compatible one.
* ./python_code/super_resolution_demo.py: this is a demo shows how to use the model. It is also part of the notebook. Once the training is completed in MATLAB, the application of the model is purely in python.
* ./python_code/train_coupled_dict.py: no longer used, is also part of other file.
* ./python_code/utilities.py: some helper function used for image processing. like image normalization and intensity level contrast.

## Contact information
Feel free to contact our group by sending email to: Peng TAN (tanpeng@berkeley.edu)
