#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#define EIGEN_DEFAULT_TO_ROW_MAJOR
#define EIGEN_DONT_PARALLELIZE
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <time.h>
#include<stdio.h>
#include<Windows.h>
#include<stdlib.h>
#include<math.h>
//#pragma comment (lib, "C:\\Program Files (x86)\\byteelSWTools\\compilers_and_libraries_2016.0.110\\windows\\compiler\\lib\\byteel64_win\\libiomp5md.lib")
using namespace std;
using namespace Eigen;

#define LINE 128
#define max 32768
#define cal_time 128
#define max_order 4   //最大迭代的阶
int  table[256];
int arc_table[256];
int inverse_table[256];
/*
配置文件格式：第1-4为切割出来的子矩阵存储的文件名，第5行是切割次数，第6行第7行是开始坐标，第8行是表示表示此次分割找到的第一个可逆
子矩阵是第几个,第9个数据是 矩阵尺寸，第10个是分隔符#
*/
struct Order                      //此结构存储迭代计算顺序
{
	int divide_time;             //切割次数
	int start_x;                //      开始坐标
	int start_y;
	int start_block;           //  表示在分割的4块中，第几块是第一个可逆子矩阵，从它开始计算
	int block_size;             //临时矩阵大小
	char block[4][128];
	//4个子矩阵的文件名
};

Order Cal_Order[cal_time];
int Read();                                     //读取配置文件函数
int Init();
int _mul(byte x, byte y);   //乘法
int _div(byte x, byte y); //除法
void Solve_1(byte *r1, byte *r2, byte mul, byte div, int size);              //size是矩阵尺寸
int GaussSolve(int n, int x, int y, byte **source, byte **result);
int ReSub(int n, int x, int y, byte  **source, byte **result);
int Cal(byte **source0, byte **source1, byte **source2, byte **source3, int rever_mark, int matrix_size, byte **the_result); //分块矩阵的方法计算矩阵的逆，前4参数为4矩阵，
																															 //第5个是可逆标记，标记第几个矩阵是第一个可逆子矩阵,matrix_ｓｉｚｅ是矩阵尺寸
int main()
{
	int start = clock();               //  计算运算时间
	char *num_buf;
	int i;
	num_buf = (char *)malloc(sizeof(char) * 4);
	FILE *source;
	int rank;
	Init();                              //初始化
	rank = Read();                     //第一步，读取配置文件
	cout << "rank: " << rank << endl;
	
	byte **temp_result[max_order];                     //存储临时结果
	byte **temp_source[4];                     //读取临时矩阵

	for (i = 0; i < rank; i++)
	{
		cout << Cal_Order[i].block[0] << endl;
		cout << Cal_Order[i].block[1] << endl;
		cout << Cal_Order[i].block[2] << endl;
		cout << Cal_Order[i].block[3] << endl;
	}

	rank = rank + 1;

	for (int d = 0; d < max_order; d++)
	{
		temp_result[d] = new byte *[max];
		for (int kc = 0; kc < max; kc++)
			temp_result[d][kc] = new byte[max];
	}                                                         //零时结果矩阵分配空间
	for (int d = 0; d < max_order; d++)
	{
		for (int kc = 0; kc < max; kc++)
		{
			for (int kr = 0; kr < max; kr++)
			{
				temp_result[d][kc][kr] = 0;
			}
		}         //临时结果数组初始化
	}
	for (i = rank - 2; i >= 0; i--)                                                 
	{
		
		cout << " I " << i << endl;
		cout << " size : " << Cal_Order[i].block_size << endl;

		for (int h = 0; h < 4; h++)
		{
			temp_source[h] = new byte *[Cal_Order[i].block_size];
			for (int f = 0; f < Cal_Order[i].block_size; f++)
			{
				temp_source[h][f] = new byte[Cal_Order[i].block_size];
			}
		}

		                                                                             //临时分配存储源矩阵的空间

		if (i == (rank - 2))
		{
			cout << "读取最初矩阵\n";
			fopen_s(&source, Cal_Order[i].block[0], "r");
			cout << Cal_Order[i].block[0] << endl;


			for (int m = 0; m < Cal_Order[i].block_size; m++)
			{
				for (int k = 0; k < Cal_Order[i].block_size; k++)
				{
					fgets(num_buf, LINE, source);
					temp_source[0][m][k] = (byte)atoi(num_buf);
				}

			}

			fclose(source);
			//读取矩阵文件
			for (int t = 0; t < Cal_Order[i].block_size; t++)
			{
				for (int h = 0; h < Cal_Order[i].block_size; h++)
				{
					temp_result[i + 1][t][h] = temp_source[0][t][h];
				}
			}

		}

		else
		{

			cout << "从上一步结果作为本步骤第一矩阵\n";
			for (int m = 0; m < Cal_Order[i].block_size; m++)
			{
				for (int k = 0; k < Cal_Order[i].block_size; k++)
				{
					temp_source[0][m][k] = temp_result[i + 1][m][k];
				}

			}

		}                                                        //如果不是最末端子矩阵，则从上一次的结果作为这一次的第一个源矩阵

		fopen_s(&source, Cal_Order[i].block[1], "r");

		for (int m = 0; m < Cal_Order[i].block_size; m++)
		{
			for (int k = 0; k < Cal_Order[i].block_size; k++)
			{
				fgets(num_buf, LINE, source);
				temp_source[1][m][k] = (byte)atoi(num_buf);
			}

		}

		fclose(source);

		fopen_s(&source, Cal_Order[i].block[2], "r");


		for (int m = 0; m < Cal_Order[i].block_size; m++)
		{
			for (int k = 0; k < Cal_Order[i].block_size; k++)
			{
				fgets(num_buf, LINE, source);
				temp_source[2][m][k] = (byte)atoi(num_buf);
			}

		}

		fclose(source);

		fopen_s(&source, Cal_Order[i].block[3], "r");


		for (int m = 0; m < Cal_Order[i].block_size; m++)
		{
			for (int k = 0; k < Cal_Order[i].block_size; k++)
			{
				fgets(num_buf, LINE, source);
				temp_source[3][m][k] = (byte)atoi(num_buf);
			}

		}

		fclose(source);
		//把4个子矩阵的数字读取到临时数组里面


		/*		for (byte t = 0; t < Cal_Order[i].block_size; t++)
		{
		for (byte h = 0; h < Cal_Order[i].block_size; h++)
		{

		temp_result[i +1][t][h] = 0;
		}
		}
		*/
		/*if (i == rank - 2)
		{
			for (int t = 0; t < Cal_Order[i].block_size; t++)
			{
				for (int h = 0; h < Cal_Order[i].block_size; h++)
				{
					if (h == t)
						temp_result[i + 1][t][h] = 1;
					else
						temp_result[i + 1][t][h] = 0;
				}
			}

			GaussSolve(Cal_Order[i].block_size, 0, 0, temp_source[0], temp_result[i + 1]);

			ReSub(Cal_Order[i].block_size, Cal_Order[i].block_size - 1, Cal_Order[i].block_size - 1, temp_source[0], temp_result[i + 1]);                                       //求单个矩阵的逆矩阵
			for (int xi = 0; xi < Cal_Order[i].block_size; xi++)
			{
				for (int yi = 0; yi < Cal_Order[i].block_size; yi++)
				{
					temp_result[i + 1][xi][yi] = _div(temp_result[i + 1][xi][yi], temp_source[0][xi][xi]);
				}
			}

		}


		*/
	
		Cal(temp_result[i + 1], temp_source[1], temp_source[2], temp_source[3], Cal_Order[i].start_block, Cal_Order[i].block_size, temp_result[i]);                 //利用公式求分块矩阵的逆
																																									//cout <<"size : " << Cal_Order[i].block_size << endl;
	
	
	
																																								/*
																																									for (byte w = 0; w<Cal_Order[i].block_size; w++)
																																									delete[] temp_result[i%2][w];                                                       //释放空间
																																									delete[] temp_result[i%2];
																																									//释放临时结果数组空间																																								*/

		for (int h = 0; h < 4; h++)
		{
			for (int w = 0; w < Cal_Order[i].block_size; w++)
				delete[] temp_source[h][w];

			delete[] temp_source[h];
		}
		//释放临时源矩阵空间

	}


	FILE *the_result;
	fopen_s(&the_result, "result.txt", "w");
	cout << "LAST  I " << i << endl;
	for (int row = 0; row < max; row++)
	{
		for (int col = 0; col < max; col++)
		{
			fprintf(the_result, "%d\n", temp_result[i + 1][row][col]);
			//cout <<" "<< temp_result[i+2][row][col]<<" ";
		}
		//cout << endl;
		//fprintf(the_result, "\n");
	}

	fclose(the_result);
	for (int d = 0; d < max_order; d++)
	{
		for (int w = 0; w < max; w++)
			delete[] temp_result[d][w];

		delete[] temp_result[d];
	}
	printf("%.3lf\n", double(clock() - start) / CLOCKS_PER_SEC);  //显示运行时间
	FILE *times;
	fopen_s(&times, "time.txt", "w");
	fprintf(times, "%.3lf\n", double(clock() - start) / CLOCKS_PER_SEC);
	fclose(times);
	cout << "over!!" << endl;
	cin.get();
	cin.get();
	cin.get();
	return 0;

}

int Init()
{
	int i = 0;
	for (i = 0; i < cal_time; i++)
	{
		Cal_Order[i].divide_time = -1;
		Cal_Order[i].start_x = -1;
		Cal_Order[i].start_y = -1;
		Cal_Order[i].start_block = -1;
	}

	table[0] = 1;//g^0  
	for (i = 1; i < 255; ++i)//生成元为x + 1  
	{
		//下面是m_table[i] = m_table[i-1] * (x + 1)的简写形式  
		table[i] = (table[i - 1] << 1) ^ table[i - 1];

		//最高指数已经到了8，需要模上m(x)  
		if (table[i] & 0x100)
		{
			table[i] ^= 0x11B;//用到了前面说到的乘法技巧  
		}
	}


	for (i = 0; i < 255; ++i)
		arc_table[table[i]] = i;


	for (i = 1; i < 256; ++i)//0没有逆元，所以从1开始  
	{
		int k = arc_table[i];
		k = 255 - k;
		k %= 255;//m_table的取值范围为 [0, 254]  
		inverse_table[i] = table[k];
	}


	return 0;
}
int Read()
{
	int len;
	int i = 0;
	char *buf;
	buf = (char *)malloc(sizeof(char) * 16);
	FILE *fp;
	fopen_s(&fp, "Config.txt", "r");
	while (!feof(fp))
	{
		fgets(buf, LINE, fp);
		len = strlen(buf);
		buf[len - 1] = '\0';
		strcpy_s(Cal_Order[i].block[0], len + 1, buf);

		fgets(buf, LINE, fp);
		len = strlen(buf);
		buf[len - 1] = '\0';
		strcpy_s(Cal_Order[i].block[1], len + 1, buf);


		fgets(buf, LINE, fp);
		len = strlen(buf);
		buf[len - 1] = '\0';
		strcpy_s(Cal_Order[i].block[2], len + 1, buf);


		fgets(buf, LINE, fp);
		len = strlen(buf);
		buf[len - 1] = '\0';
		strcpy_s(Cal_Order[i].block[3], len + 1, buf);


		fgets(buf, LINE, fp);
		Cal_Order[i].divide_time = atoi(buf);

		fgets(buf, LINE, fp);
		Cal_Order[i].start_x = atoi(buf);

		fgets(buf, LINE, fp);
		Cal_Order[i].start_y = atoi(buf);

		fgets(buf, LINE, fp);
		Cal_Order[i].start_block = atoi(buf);

		fgets(buf, LINE, fp);
		Cal_Order[i].block_size = atoi(buf);

		fgets(buf, LINE, fp);

		//过滤掉 #  分隔符
		i++;
	}
	fclose(fp);
	return i;
}
int GaussSolve(int n, int x, int y, byte **source, byte **result)  //高斯消元法求解
																   //此法编程思想：首先定义一个起点，从起点开始，把起点往下纵列元素全部改写为0，然后起点沿着对角线+1，继续迭代，直到结束
																   //并且每一步操作，都对result矩阵同样进行，直到最后返回result
{
	//cout << "test1" << endl;
	int i;
	if (n > 1 && (x + 1)<n)
	{

		//if (source[x][y] != 0)
		{
			for (i = x + 1; i < n; i++)
			{
				Solve_1(result[i], result[x], source[i][y], source[x][y], n);     //对单位矩阵进行同样操作
				Solve_1(source[i], source[x], source[i][y], source[x][y], n);    //把第二行和以后行的第一个数化为0		

			}
			return GaussSolve(n, x + 1, y + 1, source, result);                  //
																				 //迭代
		}

	}
	else
	{

		return 0;
	}

}
void Solve_1(byte *r1, byte *r2, byte mul, byte div, int size)//把r2乘以一个数加到r1上面去，为了能乘以分数，分子是mul分母是div
{
	byte temp = _div(mul, div);

	for (int y = 0; y < size; y++)                                      //下次记得修改128为矩阵尺寸
	{
		//*(r1+i) = *(r1+i) + r2[i] * mul / div;
		//*(r1 + i) = (*(r1 + i) + _div(_mul(r2[i], mul), div)) % 256;
		//*(r1 + i) = (*(r1 + i) + _mul(r2[i], temp)) % 256;
		r1[y] = (r1[y] - _mul(r2[y], temp) + 256) % 256;
	}

}
int ReSub(int n, int x, int y, byte  **source, byte **result)  //回代函数
{

	int i;


	if ((x - 1) >= 0)
	{

		for (i = x - 1; i >= 0; i--)
		{

			Solve_1(result[i], result[x], source[i][y], source[x][y], n);     //对单位矩阵进行同样操作
			Solve_1(source[i], source[x], source[i][y], source[x][y], n);    //把第二行和以后行的第一个数化为0
		}

		return ReSub(n, x - 1, y - 1, source, result);
		//迭代
	}

	else
	{
		return 0;
	}
}


int _mul(byte x, byte y)
{
	if (!x || !y)
		return 0;

	return table[(arc_table[(int)x] + arc_table[(int)y]) % 255];     //x y都不为0进行此运算
}

int  _div(byte x, byte y)
{
	return _mul((int)x, inverse_table[(int)y]);
}

int  Cal(byte **source0, byte **source1, byte **source2, byte **source3, int rever_mark, int matrix_size, byte **the_result)
{
	int i;
	int j;
	int h;
	byte ** transit = new byte *[matrix_size];                  //以下是为了求   tranist=(D-CA^(-1)B)^(-1) ;
	for (i = 0; i < matrix_size; i++)
	{
		transit[i] = new byte[matrix_size];
	}



	byte ** temp[max_order];
	for (h = 0; h< max_order; h++)
	{
		temp[h] = new byte *[matrix_size];
		for (i = 0; i < matrix_size; i++)
			temp[h][i] = new byte[matrix_size];
	}

	for (h = 0; h < max_order; h++)
		for (i = 0; i < matrix_size; i++)
			for (j = 0; j < matrix_size; j++)
				temp[h][i][j] = 0;





	//if (rever_mark == 0)                                             //如果第0个矩阵是第一个可逆子矩阵
	{

		

		//source0已经是逆矩阵，直接赋值
		//  std::array<Eigen::Matrix<double, 3, 3>, 10> aaa;
		//MatrixXf mat[12]{ 12,matrix_size,matrix_size};// = MatrixXf::Random(1000, 1000);;
		//= MatrixXf::(matrix_size);

		//	for (byte y = 0; y < 16;y++)
		//mat[y] = MatrixXf::Random(matrix_size);
		//MatrixXf Mat_A(matrix_size, matrix_size);
		// Matrix<byte,matrix_size,matrix_size > test;
		// test(1, 1) = 258;
		// cout << "test!!!" << (int)test(1, 1) << endl;
	
		//typedef Matrix<byte, Dynamic, Dynamic>  MatrixXB;
		cout << "cal start!\n";
		MatrixXi Mat_B(matrix_size,matrix_size);
		MatrixXi Mat_C(matrix_size, matrix_size);
		MatrixXi Mat_D(matrix_size, matrix_size);
		MatrixXi Mat_C_Ainv(matrix_size, matrix_size);
		MatrixXi Mat_C_Ainv_B(matrix_size, matrix_size);
		MatrixXi Mat_D_CAinvB(matrix_size, matrix_size);
		MatrixXi Mat_Tra_inv(matrix_size, matrix_size);
		MatrixXi Mat_Ainv_B(matrix_size, matrix_size);
		MatrixXi Mat_ResultA(matrix_size, matrix_size);
		MatrixXi Mat_Zero(matrix_size, matrix_size);
		MatrixXi Mat_Ainv(matrix_size, matrix_size);                    //A^(-1)
		MatrixXi Mat_ResultB(matrix_size, matrix_size);
		MatrixXi Mat_ResultC(matrix_size, matrix_size);
		MatrixXi Mat_ResultD(matrix_size, matrix_size);
																	 //Translate(source2, mat0, matrix_size);
		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{
				Mat_C(i, j) = (int)source2[i][j];
			}
		}
		//Translate(transit, mat12, matrix_size);
		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{

				Mat_Ainv(i, j) = (int)source0[i][j];
			}
		}
		//matrix_mul(source2, transit, temp[0], matrix_size);     CA^(-1)
		Mat_C_Ainv = Mat_C * Mat_Ainv;                                           //问题出在：MatrixXf 溢出了！！！！！

		//Translate(source1, mat3, matrix_size);
		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{

				Mat_B(i, j) = (int)source1[i][j];
			}
		}
		//matrix_mul(temp[0], source1, transit, matrix_size);         CA^(-1)B
		//mat3(6, 78) = (float)source1[6][78];

		Mat_C_Ainv_B = Mat_C_Ainv * Mat_B;
		//matrix_subt(source3, transit, temp[0], matrix_size);        D-CA^(-1)B
		//Translate(source3, mat4, matrix_size);
		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{

				Mat_D(i, j) = (int)source3[i][j];
			}
		}
		Mat_D_CAinvB = Mat_D - Mat_C_Ainv_B;                              //  D-CA^(-1)B

		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{
				if (i == j)
				{
					transit[i][j] = (byte)1;
					//temp[2][i][j] = 1;
					//mat1(i, j) = 1;
				//	mat5(i, j) = 1;
				}
				else
				{
					transit[i][j] = (byte)0;
					//mat1(i, j) = 0;
					//mat5(i, j) = 0;
				}
				//	temp[3] = 0;



			}

		}

		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{

				temp[0][i][j] = (byte)Mat_D_CAinvB(i, j);
			}
		}
	
		GaussSolve(matrix_size, 0, 0, temp[0], transit);                 //

		ReSub(matrix_size, matrix_size - 1, matrix_size - 1, temp[0], transit);           //公式未完善，没写完,这只是求除了第一个模块的一个公共部分，需要依照公式

		for (int xi = 0; xi < matrix_size; xi++)                                       //(D-CA^(-1)B)^(-1)
		{
			for (int yi = 0; yi < matrix_size; yi++)
			{
				transit[xi][yi] = _div(transit[xi][yi], temp[0][xi][xi]);
			}
		}                                    

		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{

				Mat_Tra_inv(i, j) = (int)transit[i][j];
			}
		}																			  //求出剩于三个模块，然后指针地址拼接起来,把4个小模块的地址指向result！！

																					  // (D-(CA ^ (-1)B)^(-1)
	
		Mat_Ainv_B = Mat_Ainv*Mat_B;
																					  //matrix_mul(temp[2], source1, temp[3], matrix_size);                    //temp3是A^(-1)B
		Mat_ResultA = Mat_Ainv + Mat_Ainv_B*Mat_Tra_inv;                              //第1模块！
		//	matrix_mul(temp[3], transit, temp[4], matrix_size);                    
	
		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{
				//the_result[i][j] = temp[4][i][j];
				the_result[i][j] = (byte)Mat_ResultA(i, j);
			}
		}

	                               	//第一个模块


	
		//matrix_subt(temp[1], temp[4], temp[4], matrix_size);           
		//Translate(temp[1], mat8, matrix_size);
		for (i = 0; i < matrix_size; i++)
		{
			for (j = 0; j < matrix_size; j++)
			{

				Mat_Zero(i, j) = 0;
			}
		}
		Mat_ResultB = Mat_Zero - Mat_Ainv_B*Mat_Tra_inv;


		for (i = 0; i <matrix_size; i++)
		{

			for (j = matrix_size; j < matrix_size * 2; j++)
			{

				//	the_result[i][j] = temp[4][i][j];                                     
				the_result[i][j] = (byte)Mat_ResultB(i , j-matrix_size);
			}
		}
		//第二个模块，拼接指针地址


		//	matrix_mul(transit, source2, temp[5], matrix_size);              //temp5=(D-CA^(-1)B)^(-1)C
		Mat_ResultC = Mat_Zero - Mat_Tra_inv*Mat_C_Ainv;
		for (i = matrix_size; i < matrix_size * 2; i++)
		{
			for (j = 0; j <matrix_size; j++)
			{
				//the_result[i][j] = temp[6][i][j];
				the_result[i][j] = (byte)Mat_ResultC(i-matrix_size, j );
			}
		}                                             //第三模块

		for (i = matrix_size; i <2 * matrix_size; i++)
		{
			for (j = matrix_size; j <2 * matrix_size; j++)
			{
				//the_result[i][j] = transit[i][j];
				the_result[i][j] = (byte)Mat_Tra_inv(i - matrix_size, j - matrix_size);
			}
		}                                //第四模块
		//cout << "mat: " << mat1(10, 13) << endl;
	}
	






	for (int d = 0; d < max_order; d++)
	{
		for (int w = 0; w < matrix_size; w++)
			delete[] temp[d][w];

		delete[] temp[d];
	}
	for (int w = 0; w <matrix_size; w++)
		delete transit[w];
	delete transit;

	return 0;
}
