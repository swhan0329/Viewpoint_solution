#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <cmath>

//setting
#define MAX 255

//#define errorlength 18588.47 //(7.7^2+7.2^2+2.8^2)^(1/3) = 118.97m -> 18588.4724pixels
#define fx 4244 //2518(1/2) //4243.807(5312X2988) 
#define fy 4229 //2513(1/2) //4229.007(5312X2988) 
#define cx 2766 //1632(1/2) //2765.728(5312X2988) 
#define cy 1510 //918(1/2) //1510.034(5312X2988) 
#define camera_ref_x 4.492 //3.9
#define camera_ref_y 1.41 //0.81 
#define camera_ref_z 1.605 //0.491 
#define feedback_max 1//5
#define D 3.0	//distance between V and O
#define PHI 3.141592

//resolution
#define RES_COL 2988. //1494.
#define RES_ROW 5312. //2656. 

//Angle of View (Android galaxy Note 5)
#define ANGLE_ROW 0.744923   //VERTICALANGLE
#define ANGLE_COL 1.2131557  //HORIZONTAL ANGLE

//view plane size
#define PLANE_ROW_SIZE 2.*D*tan(ANGLE_ROW / 2)
#define PLANE_COL_SIZE 2.*D*tan(ANGLE_COL / 2)

//#define yaw 0
//#define roll 0
//#define pitch 0

using namespace Eigen;
using namespace std;

void main()
{
	int Num_of_Space_p = 0;   // p 공백수를 저장할 변수
	int Num_of_Space_q = 0;   // q 공백수를 저장할 변수
	int statep;               // p file close
	int stateq;               // q file close
	char Stringp[MAX];        // 한 줄씩 읽어 올 문자열을 저장할 공간
	char chp;                 // p txt를 한 줄씩 읽어 온 문자열을 한 문자씩 저장
	char Stringq[MAX];        // 한 줄씩 읽어 올 문자열을 저장할 공간
	char chq;                 // q txt를 한 줄씩 읽어 온 문자열을 한 문자씩 저장

	// p의 값의 갯수를 읽고 가변적 matrix를 생성 
	FILE* filep = fopen("indoor_coordinate.txt", "rt");
	if (filep == NULL)
	{
		printf("file open error!\n");
	}
	while (1)
	{
		if (feof(filep) != 0)
			break;
		fgets(Stringp, MAX, filep);
		for (int i = 0; i < strlen(Stringp); i++)
		{
			chp = Stringp[i];
			if (chp == ' ')
				Num_of_Space_p++;
		}

		printf("공백 수 : %1d\n", Num_of_Space_p);
		statep = fclose(filep);
		if (statep != 0)
		{
			printf("file close error!\n");
			break;
		}
		break;
	}
	Num_of_Space_p += 1; // 실제 좌표는 p의 공백보다 1개 더 많음

	// p의 값의 갯수를 읽고 가변적 matrix를 생성
	FILE* fileq = fopen("pixel_coordinate.txt", "rt");
	if (fileq == NULL)
	{
		printf("file open error!\n");
	}
	while (1)
	{
		if (feof(fileq) != 0)
			break;
		fgets(Stringq, MAX, fileq);
		for (int i = 0; i < strlen(Stringq); i++)
		{
			chq = Stringq[i];
			if (chq == ' ')
				Num_of_Space_q++;
		}
		printf("공백 수 : %1d\n", Num_of_Space_q);
		stateq = fclose(fileq);
		if (stateq != 0)
		{
			printf("file close error!\n");
			break;
		}
		break;
	}
	Num_of_Space_q = Num_of_Space_q + 1; // 실제 좌표는 q의 공백보다 1개 더 많음

	//가변적 p,q 행렬 선언
	MatrixXd q = MatrixXd(Num_of_Space_q, 1); //Pixel coordinates of projected points
	MatrixXd p = MatrixXd(Num_of_Space_p, 1); //coordinates of points in 3D space
	MatrixXd R = MatrixXd((2 * Num_of_Space_p) / 3, 12);

	//p, q 행렬에 인도어 맵, 사진 좌표 데이터 넣기
	ifstream indoorFile;
	double sump = 0;
	indoorFile.open("indoor_coordinate.txt");
	if (!indoorFile.eof())
	{
		for (int i = 0; i < Num_of_Space_p; i++)
		{
			indoorFile >> p(i, 0);
		}
	}
	indoorFile.close();

	/*
	//p 입력값 확인
	cout << "p " << endl;
	for (int i = 0; i < Num_of_Space_p; i++)
	{
		cout << p(i, 0) << " ";
	}
	cout << endl;
	*/
	ifstream PixelFile;
	double sumq = 0;
	PixelFile.open("pixel_coordinate.txt");
	if (!PixelFile.eof())
	{
		for (int i = 0; i < Num_of_Space_q; i++)
		{
			PixelFile >> q(i, 0);
		}
	}
	PixelFile.close();
	/*
	//q 입력값 확인
	cout << "q " << endl;
	for (int i = 0; i < Num_of_Space_q; i++)
	{
		cout << q(i, 0) << " ";
	}
	cout << endl;
	*/
	//feedback num만큼 feedback 하기
	for (int num = 0; num < feedback_max; num++)
	{
		cout << "피드백: " << num+1 << " 번째" << endl << endl;
		//R 행렬에 데이터 집어 넣기
		for (int i = 0; i < Num_of_Space_p / 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				R(2 * i, j) = R(2 * i + 1, j + 4) = p(3 * i + j, 0);
				R(2 * i, j + 8) = (-1)*q(2 * i, 0)*p(3 * i + j, 0);
				R(2 * i + 1, j + 8) = (-1)*q(2 * i + 1, 0)*p(3 * i + j, 0);
			}
			R(2 * i, 11) = (-1)*q(2 * i, 0);
			R(2 * i + 1, 11) = (-1)*q(2 * i + 1, 0);
			R(2 * i, 4) = R(2 * i, 5) = R(2 * i, 6) = R(2 * i, 7) = 0;
			R(2 * i + 1, 0) = R(2 * i + 1, 1) = R(2 * i + 1, 2) = R(2 * i + 1, 3) = 0;
			R(2 * i, 3) = R(2 * i + 1, 7) = 1;
		}
		/*
		//R 입력값 확인
		cout << "R " << endl;
		for (int i = 0; i < 12; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				cout << R(i, j) << " ";
			}
			cout << endl;
		}
		cout << endl;
		*/

		MatrixXd RT = R.transpose(); //R 행렬을 tranpose 시킨 행렬
		/*x
		//RT 입력값 확인
		cout << "R^T " << endl;
		for (int i = 0; i < 12; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				cout << RT(i, j) << " ";
			}
			cout << endl;
		}
		*/

		MatrixXd A = RT*R; //R^T와 R의 곱행렬
		/*
		//A 입력값 확인
		cout << endl;
		cout << "R^T*R " << endl;
		for (int i = 0; i < 12; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				cout << A(i, j) << " ";
			}
			cout << endl;
		}
		cout << endl;
		*/

		EigenSolver<MatrixXd> es(A);   //A 행렬을 es로 두고
		cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl << endl; //es의 고윳값을 구한다
		cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl; //es의 고윳값을 이용하여 고유벡터를 구한다.

		float min = es.eigenvalues().coeff(0, 0).real();
		int mincoeff;
		for (int i = 0; i < 12; i++)
		{
			if (es.eigenvalues().coeff(i, 0).real() < min)
			{
				min = es.eigenvalues().coeff(i, 0).real();
				mincoeff = i;
			}
		}
		cout << endl;
		cout << "mincoeff" << mincoeff << endl; //mincoeff 값 확인
		double lambda = min; //고윳값 중 가작 작은(0에 가까운 값)을 lambda라 한다
		cout << "Consider the smallest eigenvalue, lambda = " << lambda << endl; //lambda 출력
		
		MatrixXd M = MatrixXd(3, 4); // q=pM의 M벡터
		MatrixXd MI = MatrixXd(4, 4); //p'=pM의 이론적 M벡터
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				MI(i,j)=0;
			}
		}
		
		MatrixXd Camera = MatrixXd(3, 3); //이론적 camera ex벡터
		Camera(0, 0) = fx;
		Camera(0, 1) = 0;
		Camera(0, 2) = cx;
		Camera(1, 0) = 0;
		Camera(1, 1) = fy;
		Camera(1, 2) = cy;
		Camera(2, 0) = 0;
		Camera(2, 1) = 0;
		Camera(2, 2) = 1;

		/*
		double roll2radi = 0;
		double pitch2radi = 0;
		double yaw2radi = 0;

		roll2radi = (3.141592 / 180) * roll; //x A
		pitch2radi = (3.141592 / 180) * pitch; //y B
		yaw2radi = (3.141592 / 180) * yaw; //z C

		//RxRyRz
		Ci(0, 0) = cos(yaw2radi)*cos(pitch2radi);
		Ci(0, 1) = -sin(yaw2radi)*cos(pitch2radi);
		Ci(0, 2) = sin(pitch2radi);
		Ci(0, 3) = camera_ref_x;
		Ci(1, 0) = sin(roll2radi)*sin(pitch2radi)*cos(yaw2radi) + sin(yaw2radi)*cos(roll2radi);
		Ci(1, 1) = cos(roll2radi)*cos(yaw2radi);
		Ci(1, 2) = cos(roll2radi)*sin(yaw2radi)*sin(pitch2radi) - sin(roll2radi)*cos(pitch2radi);
		Ci(1, 3) = camera_ref_y;
		Ci(2, 0) = sin(roll2radi)*sin(yaw2radi)*cos(pitch2radi)-cos(roll2radi)*sin(pitch2radi);
		Ci(2, 1) = sin(roll2radi)*cos(yaw2radi);
		Ci(2, 2) = sin(roll2radi)*sin(yaw2radi)*sin(pitch2radi)+cos(roll2radi)*cos(pitch2radi);
		Ci(2, 3) = camera_ref_z;

		//RyRxRz
		Ci(0, 0) = cos(yaw2radi)*cos(pitch2radi) - sin(yaw2radi)*sin(roll2radi)*sin(pitch2radi);
		Ci(0, 1) = -sin(yaw2radi)*cos(roll2radi);
		Ci(0, 2) = cos(yaw2radi)*sin(pitch2radi) + sin(yaw2radi)*sin(roll2radi)*cos(pitch2radi);
		Ci(0, 3) = camera_ref_x;
		Ci(1, 0) = sin(yaw2radi)*cos(pitch2radi) + cos(yaw2radi)*sin(roll2radi)*sin(pitch2radi);
		Ci(1, 1) = cos(yaw2radi)*cos(roll2radi);
		Ci(1, 2) = sin(yaw2radi)*sin(pitch2radi) - cos(yaw2radi)*sin(roll2radi)*cos(pitch2radi);
		Ci(1, 3) = camera_ref_y;
		Ci(2, 0) = -cos(roll2radi)*sin(pitch2radi);
		Ci(2, 1) = sin(roll2radi);
		Ci(2, 2) = cos(roll2radi)*cos(pitch2radi);
		Ci(2, 3) = camera_ref_z;


		cout << "이론 R|T" << endl << endl;
		cout << Ci << endl;
		cout << "이론 intrinsic matrix" << endl << endl;
		cout << Camera << endl<<endl;
		MI = 0.0001*Camera*Ci;
		cout << "이론 M" << endl << MI<<endl;
		*/
		/*
		//lambda에 해당하는 고유벡터 출력 확인
		for (int i = 0; i < 12; i++)
		{
		for (int j = 0; j < 12; j++)
		{
		cout << es.eigenvectors().coeff(i, j).real() << "   ";
		}
		cout << endl;
		}
		cout << endl;
		*/

		//M벡터에 고유 벡터 coefficient 입력
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				M(i, j) = es.eigenvectors().coeff(4 * i + j, mincoeff).real();
			}
		}

		// cout << "Here is the matrix 계산 M:" << endl << M << endl << endl;

		MatrixXd V = M.fullPivLu().kernel(); //MV=0을 풀기 위해 null space of M을 구한다.
		cout << "Here is a solution V to the equation MV=0:" << endl << endl << V << endl << endl; //homogenous coordinate V 좌표를 구한다.
		cout << "Here is a Viewpoint coordinate(x, y, z): " << "( " << V(0, 0) / V(3, 0) << ", " << V(1, 0) / V(3, 0) << ", " << V(2, 0) / V(3, 0) << " )" << endl; //실제 3-D좌표계로 변환
		
		//feedback 추가 2017.12.27
		MatrixXd n = MatrixXd(1, 4);
		MatrixXd tempR= MatrixXd(1, 3);
		MatrixXd tempS = MatrixXd(1, 3);
		MatrixXd n_vec = MatrixXd(1, 3);  // view plane normal vector
		MatrixXd R = MatrixXd(1, 4); // view plane X-axes
		MatrixXd S = MatrixXd(1, 4); // view plane Y-axes
		MatrixXd O = MatrixXd(1, 4); // view plane Origin
		
		double rot_angle = 0.0; //rotated angle(degree)
		n_vec(0, 0) = 1.;  n_vec(0, 1) = 0.;  n_vec(0, 2) = 0.;
		double mother = sqrt(pow(n_vec(0, 0), 2) + pow(n_vec(0, 1), 2) + pow(n_vec(0, 2), 2));

		for (int i = 0; i < 3; i++)
		{
			//normalize n_vector
			n_vec(0, i) /= mother;
			//set the veiw plane equation n
			n(0, i) = n_vec(0, i);
			//set the veiw plane Origin
			O(0, i) = V(i, 0) / V(3, 0) + (D*n_vec(0, i));
		}

		/////////////////////////No Rotation//////////////////////////////
		O(0, 3) = 1;
		n(0, 3) = -1 * ((n(0, 0) * V(0, 0) / V(3, 0)) + (n(0, 1) * V(1, 0) / V(3, 0)) + (n(0, 2) * V(2, 0) / V(3, 0))) - D;
		//printf("\nn = {%5.5f, %5.5f, %5.5f, %5.5f}\n", n(0,0), n(0,1), n(0,2), n(0,3));
		//printf("O = {%10.2f, %10.2f, %10.2f, %10.2f}\n", O(0), O(1), O(2), O(3));

		//R Setting
		R(0, 0) = n(0, 1);
		R(0, 1) = (-1)*n(0, 0);
		R(0, 2) = 0.;

		double temp = sqrt(pow(n(0, 0), 2) + pow(n(0, 1), 2));
		double s_angle = (PHI / 2) - abs(asin(n(0, 2)));
		double temp_sz = temp*tan(s_angle);
		
		cout << "temp" << temp << endl;
		cout << "s_angle" << s_angle << endl;
		cout << "R" << endl << R << endl;
		cout << "O" << endl << O << endl;
		
		//S vector setting
		if (n(0,2) == 0.)
		{
			S(0,0) = 0.;
			S(0,1) = 0.;
			S(0,2) = 1.;
		}
		else if (n(0,2) > 0.)
		{
			S(0,0) = -1 * n(0,0);
			S(0,1) = -1 * n(0,1);
			S(0,2) = temp_sz;
		}
		else
		{
			S(0,0) = n(0,0);
			S(0,1) = n(0,1);
			S(0,2) = temp_sz;
		}
		
		//////////////////////////////////////////Normalized////////////////////////////////////////////////
		double mom1 = sqrt(pow(R(0,0), 2) + pow(R(0,1), 2) + pow(R(0,2), 2));
		double mom2 = sqrt(pow(S(0,0), 2) + pow(S(0,1), 2) + pow(S(0,2), 2));
		
		for (int i = 0; i < 3; i++)
		{
			R(0,i) /= mom1;
			S(0,i) /= mom2;
		}
		R(0,3) = 0.;
		S(0,3) = 0.;
		//printf("R0 = {%10.2f, %10.2f, %10.2f, %10.2f}\nS = {%10.2f, %10.2f, %10.2f, %10.2f}\n\n", R(0), R(1), R(2), R(3), S(0), S(1), S(2), S(3));

		//////////////////////////////////////////Rotation//////////////////////////////////////////////////
		for (int i = 0; i < 3; i++)
		{
			tempR(0,i) = R(0,i);
			tempS(0,i) = S(0,i);
		}
		printf("rot angle=%f\n", rot_angle);
		double rotation = rot_angle*PHI / 180.;	//rotation angle setting
		
	    //rotation adaptation
		for (int i = 0; i < 3; i++)
		{
			if ((rotation > 0 && rotation < PHI / 2.) || (rotation > 3. * PHI / 2. && rotation < 2. * PHI))
			{
				R(0,i) = tempR(0,i) + (tempS(0,i) * tan(rotation));
				S(0,i) = tempS(0,i) - (tempR(0,i) * tan(rotation));
			}
			else if (rotation > PHI / 2.&&rotation < 3. * PHI / 2.)
			{
				R(0,i) = (-1.*tempR(0,i)) - (tempS(0,i) * tan(rotation));
				S(0,i) = (-1.*tempS(0,i)) + (tempR(0,i) * tan(rotation));
			}
			else if (rotation == PHI / 2.)
			{
				R(0,i) = tempS(0,i);
				S(0,i) = -1.*tempR(0,i);
			}
			else if (rotation == PHI)
			{
				R(0,i) = -1 * tempR(0,i);
				S(0,i) = -1.*tempS(0,i);
			}
			else if (rotation == 3.*PHI / 2.)
			{
				R(0,i) = -1 * tempS(0,i);
				S(0,i) = tempR(0,i);
			}
		}

		//normalized
		mom1 = sqrt(pow(R(0,0), 2) + pow(R(0,1), 2) + pow(R(0,2), 2));
		mom2 = sqrt(pow(S(0,0), 2) + pow(S(0,1), 2) + pow(S(0,2), 2));

		for (int i = 0; i < 3; i++)
		{
			R(0,i) /= mom1;
			S(0,i) /= mom2;
		}
		R(0, 3) = 0.;
		S(0, 3) = 0.;
		cout << "R" << endl << R << endl << endl;
		cout << "S" << endl << S << endl << endl;

		MI(0, 0) = (-1)*n(0, 1)*V(1, 0)/V(3, 0) - n(0, 2)* V(2, 0) / V(3, 0) - n(0, 3)*V(3, 0/V(3, 0)) / V(3, 0);
		MI(0, 1) = n(0, 0)*V(1, 0) / V(3, 0);
		MI(0, 2) = n(0, 0)*V(2, 0) / V(3, 0);
		MI(0, 3) = n(0, 0)*V(3, 0) / V(3, 0);
		MI(1, 0) = n(0, 1)*V(0, 0) / V(3, 0);
		MI(1, 1) = (-1)*n(0, 0)*V(0, 0) / V(3, 0) - n(0, 2)*V(2, 0) / V(3, 0) - n(0, 3)*V(3, 0) / V(3, 0);
		MI(1, 2) = n(0, 1)*V(2, 0) / V(3, 0);
		MI(1, 3) = n(0, 1)*V(3, 0) / V(3, 0);
		MI(2, 0) = n(0, 2)*V(0, 0) / V(3, 0);
		MI(2, 1) = n(0, 2)*V(1, 0) / V(3, 0);
		MI(2, 2) = (-1)*n(0, 0)*V(0, 0) / V(3, 0) - n(0, 1)*V(1, 0) / V(3, 0) - n(0, 3)*V(3, 0) / V(3, 0);
		MI(2, 3) = n(0, 2)*V(3, 0) / V(3, 0);
		MI(3, 0) = n(0, 3)*V(0, 0) / V(3, 0);
		MI(3, 1) = n(0, 3)*V(1, 0) / V(3, 0);
		MI(3, 2) = n(0, 3)*V(2, 0) / V(3, 0);
		MI(3, 3) = (-1)*n(0, 0)*V(0, 0) / V(3, 0) - n(0, 1)*V(1, 0) / V(3, 0) - n(0, 2)*V(2, 0) / V(3, 0);
		
		//MI 헁렬을 이용하여 모든 p점을 q로 변환한다.
		MatrixXd newp = MatrixXd(6, 4);
		newp(0, 0)= p(0,0); newp(0, 1)= p(1, 0); newp(0, 2)= p(2, 0); newp(0, 3)=1;
		newp(1, 0)= p(3, 0); newp(1, 1)= p(4, 0); newp(1, 2)= p(5, 0); newp(1, 3)=1;
		newp(2, 0)= p(6, 0); newp(2, 1)= p(7, 0); newp(2, 2)= p(8, 0); newp(2, 3)=1;
		newp(3, 0)= p(9, 0); newp(3, 1)= p(10, 0); newp(3, 2)= p(11, 0); newp(3, 3)=1;
		newp(4, 0)= p(12, 0); newp(4, 1)= p(13, 0); newp(4, 2)= p(14, 0); newp(4, 3)=1;
		newp(5, 0)= p(15, 0); newp(5, 1)= p(16, 0); newp(5, 2)= p(17, 0); newp(5, 3)=1;
		
		MatrixXd p_1 = MatrixXd(6, 4); //p'
		//cout << MI << endl;
		p_1 = newp * MI;
		//cout << Compare; 
	
		MatrixXd p_2 = MatrixXd(6, 3); //p''=p'VC
		MatrixXd K = MatrixXd(3, 4);
		MatrixXd VC = MatrixXd(4, 3);
		K(0, 0) = R(0, 0); K(0, 1) = R(0, 1); K(0, 2) = R(0, 2); K(0, 3) = 0;
		K(1, 0) = S(0, 0); K(1, 1) = S(0, 1); K(1, 2) = S(0, 2); K(1, 3) = 0;
		K(2, 0) = O(0, 0); K(2, 1) = O(0, 1); K(2, 2) = O(0, 2); K(2, 3) = 1;
		
		cout << "K" << endl << K << endl << endl;
		VC = K.transpose()*(K*K.transpose()).inverse();
		cout << "VC" << endl << VC << endl << endl;
		p_2 = p_1*VC;
	
		MatrixXd p_3 = MatrixXd(6, 3); //P'''=P''VP
		double temp_plane_row = PLANE_ROW_SIZE;
		double temp_plane_col = PLANE_COL_SIZE;
		double var_row = (double)RES_COL / temp_plane_row;
		double var_col = (double)RES_ROW / temp_plane_col;

		for (int i = 0; i < 6; i++) 
		{
			p_3(i, 0) = p_2(i, 0)*var_col;
			p_3(i, 1) = p_2(i, 1)*var_row;
			p_3(i, 2) = 1;
		}

		cout << endl;
		cout <<" p_1" << endl << p_1 << endl << endl;
		cout <<" p_2" << endl << p_2 << endl << endl;
		cout <<" p_3 " << endl << p_3 << endl << endl;
		MatrixXd pixels = MatrixXd(6, 3);
		for (int i = 0; i < 6; i++)
		{
			pixels(i, 0) =( p_3(i, 0)*((int)RES_COL / 2));
			pixels(i, 1) = ((p_3(i, 1))*(-1) + ((int)RES_ROW / 2));
			pixels(i, 2) = 1;
		}
		cout << "pixels " << endl << pixels << endl << endl;
		//error function
		double error_sum = 0;
	}
}