/*=========================================================================

Program:   SIFT3D : sift 3D with vtk readers
Module:    SIFT3D
Language:  C++
Date:      2017/06
Auteur:   Sebastien Valette

=========================================================================*/

#include <sstream>
#include <vtkImageCast.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNIFTIImageWriter.h>

#include "vtkRobustImageReader.h"

#include "immacros.h"
#include "imutil.h"
#include "sift.h"

/* Message buffer size */
#define BUF_SIZE 1024


void err_msgu( const char * text ) {
	std::cout<<text<<std::endl;
}

#define HIST_GET_IDX(a, p) ((a) + (p) * NBINS_AZ)
#define HIST_GET(hist, a, p) ((hist)->bins[HIST_GET_IDX(a, p)])

int main( int argc, char *argv[] )
{

	if (argc<2)
	{
		cout<<"Usage : SIFT3D volume -s spacing -ct cornerThreshold -pt peakThreshold"<<endl;
		cout<<"options : "<<endl;
		cout<<"peakThreshold  : values between 0 and 1, default value = 0.1"<<endl;
		cout<<"cornerThreshold  : values between 0 and 1, default value = 0.5"<<endl;
		cout<<"spacing  : resample input volume with isotropic spacing"<<endl;
		exit(1);
	}

	SIFT3D sift3d;
	double spacing[ 3 ];
	Image imReal;
	Image *im = &imReal;
	Keypoint_store kp;
	SIFT3D_Descriptor_store desc;

	// Initialize the SIFT data 
	if (init_SIFT3D(&sift3d)) {
		std::cout <<  "Failed to initialize SIFT data." << std::endl;
		exit( 1 );
	}

	// Initialize data 
	init_Keypoint_store( &kp ); 
	init_SIFT3D_Descriptor_store( &desc );
	init_im( im );

	// Parse optionnal arguments
	int argumentsIndex = 2;
	double desiredSpacing = -1;
	int writeCSV = 0;
	int drawPoints = 0;

	while ( argumentsIndex < argc ) {

		char* key = argv[ argumentsIndex ];
		char *value = argv[ argumentsIndex + 1 ];

		if (strcmp(key, "-s") == 0) {

			desiredSpacing = atof( value );
			std::cout << "desired spacing : " << desiredSpacing << std::endl;

		}

		if (strcmp(key, "-pt") == 0) {

			std::cout << "peak threshold : " << value << std::endl;
			sift3d.peak_thresh = atof( value );

		}

		if (strcmp(key, "-ct") == 0) {

			std::cout << "corner threshold : " << value << std::endl;
			sift3d.corner_thresh = atof( value );

		}


		if (strcmp(key, "-dp") == 0) {

			std::cout << "draw points : " << value << std::endl;
			drawPoints = atoi( value );

		}

		if (strcmp(key, "-csv") == 0) {

			std::cout << "write csv files : " << value << std::endl;
			writeCSV = atoi( value );

		}

		argumentsIndex += 2;
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	vtkRobustImageReader *Reader = vtkRobustImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();

	int dimensions[3];
	Reader->GetOutput()->GetDimensions(dimensions);
	cout << "dmensions : " << dimensions[ 0 ] << " " << dimensions[1] 
		<< " " << dimensions[ 2 ] << std::endl;

	vtkImageCast *cast = vtkImageCast::New();
	cast->SetInputData( Reader->GetOutput() );
	cast->SetOutputScalarType( VTK_FLOAT );
	cast->Update();

	vtkImageData *resampled;

	if (desiredSpacing > 0 ) {

		std::cout << "new spacing : " << desiredSpacing << std::endl;

		vtkImageResample *Resampler=vtkImageResample::New();
		Resampler->SetInputData( cast->GetOutput() );
		double *sp = Reader->GetOutput()->GetSpacing();

		for (int i = 0; i < 3; i++) {

			Resampler->SetAxisMagnificationFactor( i, sp[ i ] / desiredSpacing );

		}

		std::cout << "Resampling...";
		Resampler->Update();
		cout<<"...done" << std::endl;
		resampled = Resampler->GetOutput();
		vtkNIFTIImageWriter *writer = vtkNIFTIImageWriter::New();
		writer->SetInputData( resampled );
		writer->SetFileName( "resampled.nii.gz" );
		writer->Write();

		resampled->GetDimensions(dimensions);
		cout << "new dmensions : " << dimensions[ 0 ] << " " << dimensions[1] 
			<< " " << dimensions[ 2 ] << std::endl;
	} else {

		resampled = cast->GetOutput();
		
	}

	resampled->GetDimensions( dimensions );
	resampled->GetSpacing( spacing );
	double scaleSpacing = pow( spacing[ 0 ] * spacing[ 1 ] * spacing[ 2 ], 1.0 / 3.0 );


	// Store the real world coordinates
	im->ux = spacing[ 0 ];
	im->uy = spacing[ 1 ];
	im->uz = spacing[ 2 ];

	// Resize im    
	im->nx = dimensions[ 0 ];
	im->ny = dimensions[ 1 ];
	im->nz = dimensions[ 2 ];
	im->nc = 1;
	im_default_stride(im);
	im_resize(im);

	// copy data
	float *in = ( float* ) resampled->GetScalarPointer();
	float *out = im->data;
	unsigned int nVoxels = dimensions[ 0 ]* dimensions[ 1 ] * dimensions[ 2 ];

	for ( unsigned int i = 0 ; i < nVoxels; i++ ) {

		out[ i ] = in[ i ];

	}

	im_scale( im );

	// Extract keypoints
	if (SIFT3D_detect_keypoints(&sift3d, im, &kp)) {

		err_msgu("Failed to detect keypoints.");
		return 1;

	}

	// Extract descriptors
	if (SIFT3D_extract_descriptors(&sift3d, &kp,&desc)) {

		err_msgu("Failed to extract descriptors.");
		return 1;

	}

	if ( writeCSV ) {

		// Optionally write the keypoints 
		if (write_Keypoint_store( "keys.csv", &kp)) {

				err_msgu("Failed to write the keypoints ");
				return 1;
		}

		// Write the descriptors
		if (write_SIFT3D_Descriptor_store("desc.csv", &desc)) {

				err_msgu("Failed to write the descriptors");
				return 1;

		}

	}

	if ( drawPoints ) {

		// Optionally draw the keypoints
		Image draw;
		Mat_rm keys;

		// Initialize intermediates
		init_im(&draw);

		if (init_Mat_rm(&keys, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE))
				err_msgu("Failed to initialize keys matrix");

		// Convert to matrices
		if (Keypoint_store_to_Mat_rm(&kp, &keys)) {
				err_msgu("Failed to convert the keypoints to "
						 "a matrix.");
				return 1;
		}

		// Draw the points
		if (draw_points(&keys, SIFT3D_IM_GET_DIMS(im), 1, &draw)) {
				err_msgu("Failed to draw the points.");
				return 1;
		}

		draw.ux = spacing[ 0 ];
		draw.uy = spacing[ 1 ];
		draw.uz = spacing[ 2 ];

		// copy data
		out= ( float* ) resampled->GetScalarPointer();
		in = draw.data;

		for ( unsigned int i = 0 ; i < nVoxels; i++ ) {
			out[i] = in[i];
		}
		vtkNIFTIImageWriter *writer = vtkNIFTIImageWriter::New();
		writer->SetInputData( resampled );
		writer->SetFileName( "points.nii.gz" );
		writer->Write();
		// Clean up
		im_free(&draw);
	}


	float valF;
	char valC;

	ofstream pointsFile;
	pointsFile.open( "points.csv", std::ofstream::out | std::ofstream::trunc);

    Mat_rm mat;
	int i, i_R, j_R;

	const int num_rows = desc.num;

	// Write the keypoints 

	double *origin = resampled->GetOrigin();

	std::cout << num_rows << " points" << std::endl;

	// Initialize the matrix
	init_Mat_rm(&mat, 0, 0, SIFT3D_FLOAT, SIFT3D_FALSE);

	// Write the data into the matrix 
	SIFT3D_Descriptor_store_to_Mat_rm(&desc, &mat);


	for (i = 0; i < num_rows; i++) {
//	std::cout << i << " point" << std::endl;

		const SIFT3D_Descriptor *const des = desc.buf + i;

		// Write the coordinates 
		pointsFile << des->xd * spacing[ 0 ] + origin[0] << ',';

		pointsFile << des->yd * spacing[ 1 ] + origin[1] << ',';

		pointsFile << des->zd * spacing[ 2 ] + origin[2] << ',';

		pointsFile << std::fixed << des->sd * scaleSpacing << ',';

		pointsFile << "1,"; //point.laplacian;

		pointsFile << std::fixed << des->rd; //point.response;

		// write the feature vector
		for ( unsigned int j = 0; j < 768 ; j++) {

			pointsFile << std::fixed << SIFT3D_MAT_RM_GET( &mat, i, j + 3, float );
			if ( j <767 ) pointsFile << ",";

		}

		pointsFile << endl;

	}

	pointsFile.close();

}
