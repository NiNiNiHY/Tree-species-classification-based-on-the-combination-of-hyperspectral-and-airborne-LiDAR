#include <liblas/liblas.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <boost/thread/thread.hpp>

#include <float.h>
/*#include <mat.h> */

#include <fstream>    //�ļ�����ͷ�ļ�

struct PointCloud
{
	double x;
	double y;
	double z;
	int intensity;
	int r;
	int g;
	int b;
	int ReturnNumber;
	int NumberOfReturn;
	int tree = 0;
	double Canopy_ratio = 0;
	double crown = 0;
}AllPoint[653491], Point[395643], PointRedundancy[36400];
//AllPoint ���е���  Point�״λز�����  

struct TreePoint
{
	double x;
	double y;
	double z;
	int X = 0;//�߽��ǩ
	int tree = 0; //����ǩ
	int maxflag = 0; //������ֵ
	int localmaxpoint = 1; //�ֲ������ǩ
	double Canopy_ratio = 0;
	double crown = 0;  //�ڷ�
} TreePoint[36400], TreePoint1[36400];//ȥ�����

double lidarmat[40][40][44] = { 0 };
//int matsize = 101;


int Allcount = 0, count = 0, number = 0;
double a = 0.0;
int effective = 0;  //��Ч�����
int maxx, minx, maxy, miny, maxz, minz;
int countR = 0;
int X = 0;

double MultipleY = 0;
double MultipleX = 0;

void MatrixTransform(double number, int tree, int Z, double feature)
{
	int X, Y;
	int i;
	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{
			X = AllPoint[i].x * MultipleX;
			Y = AllPoint[i].y * MultipleY;
			lidarmat[X][39 - Y][Z] = feature;
		}
	}
}
void H_Matrix(double number)  //��������ƽ���߶�
{
	int X, Y;
	int i, j;
	int tree_count[40][40] = { 0 };
	int tree_count1[40][40] = { 1 };

	for (i = 0; i < number; i++)
	{

		X = AllPoint[i].x * MultipleX;
		Y = AllPoint[i].y * MultipleY;
		if(AllPoint[i].b!=0) lidarmat[X][39 - Y][0] += 1 ;
		lidarmat[X][39 - Y][1] += AllPoint[i].z / 100.0;
		lidarmat[X][39 - Y][2] += AllPoint[i].intensity / 100.0;
		if (AllPoint[i].ReturnNumber == 1)
		{
			lidarmat[X][39 - Y][2] += AllPoint[i].intensity / 100.0;
			tree_count1[X][39 - Y]++;
		}
		tree_count[X][39 - Y]++;

	}
	for (i = 0; i < 40; i++)
		for (j = 0; j < 40; j++)
		{
			
			lidarmat[i][j][1] = (lidarmat[i][j][1] / tree_count[i][j])*100.0;
			lidarmat[i][j][2] = (lidarmat[i][j][2])*100.0;  //Isum
			lidarmat[i][j][3] = (lidarmat[i][j][2] / tree_count[i][j])*100.0; //Imean
			//lidarmat[i][j][6] = (lidarmat[i][j][1])*100.0;  //Isum
			//lidarmat[i][j][5] = (lidarmat[i][j][1] / tree_count[i][j])*100.0; //Imean
		}
}


double max(char xyz, int count)

{
	//max�����β�Ϊ�ṹ�����飬ʵ��Ϊ�ṹ��ָ�롣 
	int  i, m = 0;
	double pointmax = 0;
	for (i = 0; i<count; i++)
	{
		switch (xyz)
		{
		case 1:if (Point[i].x>pointmax) pointmax = Point[i].x; break;
		case 2:if (Point[i].y>pointmax) pointmax = Point[i].y; break;
		case 3:if (Point[i].z>pointmax) pointmax = Point[i].z; break;
		}
	}
	return pointmax;
}
double min(char xyz, int count)

{
	//max�����β�Ϊ�ṹ�����飬ʵ��Ϊ�ṹ��ָ�롣 
	int  i, m = 0;
	double pointmin = 10000000;
	for (i = 0; i<count; i++)
	{
		switch (xyz)
		{
		case 1:if (Point[i].x<pointmin) pointmin = Point[i].x; break;
		case 2:if (Point[i].y<pointmin) pointmin = Point[i].y; break;
		case 3:if (Point[i].z<pointmin) pointmin = Point[i].z; break;
		}
	}
	return pointmin;
}

int EliminationOfRedundancy(double a)  //ȥ����
{
	double xl = 0, yl = 0;
	double xh = a, yh = a;
	int maxpointz = 0;
	int i = 0;
	countR = 0;
	int countr = 0;

	for (xl = 0; xl < maxx; xl += a)
	{
		for (yl = 0; yl < maxy; yl += a)
		{
			countr++;
		}
	}

	countr = 0;
	for (xl = 0; xl < maxx; xl += a)
		for (yl = 0; yl < maxy; yl += a)
		{
			for (i = 0; i < count; i++)
			{
				if (((Point[i].x< (xl + a)) && (Point[i].x >= xl)) && ((Point[i].y<(yl + a)) && (Point[i].y >= yl)))
				{
					if (maxpointz < Point[i].z)
					{
						maxpointz = Point[i].z;
						PointRedundancy[countR] = Point[i];
						TreePoint[countR].x = xl;
						TreePoint[countR].y = yl;
						TreePoint[countR].z = Point[i].z;
						TreePoint1[countR] = TreePoint[countR];
					}
				}
			}
			maxpointz = 0;
			countR++;
		}
	cout << "treepoint=" << countR << endl;
	return countR;
}

void Filter(double a)  //�˲�
{
	int i, j, k, p;
	double core[3][3] = { 0 };
	double sum;
	int flag = 0;
	for (i = 0; i < countR; i++)  //��ȡ��ǰ�����Ҹ���Ĵ�С
	{
		for (k = 0; k < 3; k++)
			for (p = 0; p < 3; p++)
			{
				core[k][p] = 0;
			}
		sum = 0;
		for (j = 0; j < countR; j++)  //��ȡ��ǰ�����Ҹ���Ĵ�С
		{
			if ((TreePoint[j].x == TreePoint[i].x) && (TreePoint[j].y == TreePoint[i].y)) { core[1][1] = 3.0* TreePoint[j].z; }//��
			if ((TreePoint[j].x == (TreePoint[i].x - a)) && (TreePoint[j].y == TreePoint[i].y)) { core[1][0] = 0.5* TreePoint[j].z; }//��
			if ((TreePoint[j].x == (TreePoint[i].x + a)) && (TreePoint[j].y == TreePoint[i].y)) { core[1][2] = 0.5* TreePoint[j].z; }//��
			if ((TreePoint[j].x == TreePoint[i].x) && (TreePoint[j].y == (TreePoint[i].y + a))) { core[0][1] = 0.5 *TreePoint[j].z; }//��
			if ((TreePoint[j].x == TreePoint[i].x) && (TreePoint[j].y == (TreePoint[i].y - a))) { core[2][1] = 0.5* TreePoint[j].z; }//��
			if ((TreePoint[j].x == (TreePoint[i].x - a)) && (TreePoint[j].y == (TreePoint[i].y + a))) { core[0][0] = 0.25*TreePoint[j].z; }//����
			if ((TreePoint[j].x == (TreePoint[i].x - a)) && (TreePoint[j].y == (TreePoint[i].y - a))) { core[2][0] = 0.25*TreePoint[j].z; }//����
			if ((TreePoint[j].x == (TreePoint[i].x + a)) && (TreePoint[j].y == (TreePoint[i].y + a))) { core[0][2] = 0.25* TreePoint[j].z; }//����
			if ((TreePoint[j].x == (TreePoint[i].x + a)) && (TreePoint[j].y == (TreePoint[i].y - a))) { core[2][2] = 0.25*TreePoint[j].z; }//����
		}

		for (k = 0; k < 3; k++)
			for (p = 0; p < 3; p++)
			{
				//cout << core[k][p] << endl;
				sum += core[k][p];
				//if (core[k][p] == 0) flag = 1;  //����ٽ�ֵȱʧ������ԭֵ
			}
		//cout << flag << endl;
		TreePoint[i].z = sum / 6.0;
		TreePoint1[i].z = TreePoint[i].z;
		//flag = 0;
	}

}

void LocalMaxPoint()  //Ѱ�Ҿֲ����ֵ ,���ؾֲ�����ֵ��ĸ���
{
	int  i = 0, j = 0, p = 0, k = 0;
	int LocalPointNumber = 0;
	double core[3][3] = { 0 };
	for (i = 0; i < countR; i++)
	{
		for (j = 0; j < countR; j++)
		{
			if ((TreePoint[j].x == TreePoint[i].x) && (TreePoint[j].y == TreePoint[i].y)) { core[1][1] = TreePoint[j].z; }//��
			if ((TreePoint[j].x == (TreePoint[i].x - a)) && (TreePoint[j].y == TreePoint[i].y)) { core[1][0] = TreePoint[j].z; }//��
			if ((TreePoint[j].x == (TreePoint[i].x + a)) && (TreePoint[j].y == TreePoint[i].y)) { core[1][2] = TreePoint[j].z; }//��
			if ((TreePoint[j].x == TreePoint[i].x) && (TreePoint[j].y == (TreePoint[i].y + a))) { core[0][1] = TreePoint[j].z; }//��
			if ((TreePoint[j].x == TreePoint[i].x) && (TreePoint[j].y == (TreePoint[i].y - a))) { core[2][1] = TreePoint[j].z; }//��
			if ((TreePoint[j].x == (TreePoint[i].x - a)) && (TreePoint[j].y == (TreePoint[i].y + a))) { core[0][0] = TreePoint[j].z; }//����
			if ((TreePoint[j].x == (TreePoint[i].x - a)) && (TreePoint[j].y == (TreePoint[i].y - a))) { core[2][0] = TreePoint[j].z; }//����
			if ((TreePoint[j].x == (TreePoint[i].x + a)) && (TreePoint[j].y == (TreePoint[i].y + a))) { core[0][2] = TreePoint[j].z; }//����
			if ((TreePoint[j].x == (TreePoint[i].x + a)) && (TreePoint[j].y == (TreePoint[i].y - a))) { core[2][2] = TreePoint[j].z; }//����
		}
		for (p = 0; p < 3; p++)
			for (k = 0; k < 3; k++)
			{
				if (core[1][1] < core[p][k]) { TreePoint[i].localmaxpoint = 0; }

			}
	}

	for (i = 0; i < countR; i++)
	{
		if (TreePoint[i].localmaxpoint == 1) { cout << i << endl; LocalPointNumber++; }
	}

}
/********************************
�����б�����
������ TreeZ ��������С���Ķ���߶ȣ�С�ڴ˸߶�ֵ������Ϊ��һ����
model ��ģʽѡ��
model = 1 ��ʾ��Ѱ�Ҿֲ�����ֵ���ж�Ϊ�������ٽ��������½���Ѱ���ڡ�
model = 0 ��ʾ�����ֵ��ʼ�����������½�����Ѱ������Ȼ��Ѱ��δʶ������ֵ�������������½���

********************************/

int TrendOfDiscriminant(double TreeZ, bool model)  //�����б�  ����������Сֵ����ֵ
{
	int treenumber = 0;
	int i = 0, j = 0, k = 0, p = 0;
	double diameterlf, diameterud, diametermax;  //����ֱ��������ֱ��
	bool flagmax = 1;//���ֵ��Ѱ������ֵ
					 //for (treenumber = 1; treenumber < TreeNumber+1; treenumber++)
					 //{
	if (model == 1) LocalMaxPoint(); //Ѱ�Ҿֲ����ֵ

	while (flagmax == 1)
	{
		treenumber++;
		i = 0, j = 0, k = 0, p = 0;
		double z = 0;
		int numbermax = 0;  //���ֵ��
		int processpoint[3] = { 0,0,0 };
		double left, right, up, down;
		bool flag = 0;
		bool flaglr = 0;  //������Ѱ�ɹ���־
		int direction = 0;// 0���� 1����  2���� 3����
		int treelift, treeright;
		int treeup[100], treedown[100];
		if (model == 0)
		{
			for (i = 0; i < countR; i++)
			{
				if ((TreePoint[i].z > z) && (TreePoint[i].tree == 0))
				{
					z = TreePoint[i].z;
					numbermax = i;
					flagmax = 1; //���ֵ�����ɹ�
				}
			}
			if (z <= TreeZ) flagmax = 0;  //�����⵽�����ֵС�ڵ���5 ���������
		}

		if (model == 1)
		{
			flagmax = 0;
			for (i = 0; i < countR; i++)
			{
				if ((TreePoint[i].localmaxpoint == 1) && (TreePoint[i].tree == 0) && (TreePoint[i].z >= TreeZ))
				{
					z = TreePoint[i].z;
					numbermax = i;
					flagmax = 1; //���ֵ�����ɹ�

					break;
				}
			}
			if (flagmax == 0) break;
		}

		//if (z <= 10) break;  //�����⵽�����ֵС�ڵ���5 ���������
		//
		//cout << "��ʼ��XΪ��" << endl;
		processpoint[1] = numbermax;
		cout << "���ֵ�㣺" << processpoint[1] << endl;
		cout << "���ֵ�߶ȣ�" << TreePoint[processpoint[1]].z << endl;
		while (direction < 2)  //x���ұ߽��
		{
			//if (flag == 1) { processpoint = numbermax; flag = 0; }
			for (i = 0; i < countR; i++)  //��ȡ��ǰ�����Ҹ���Ĵ�С
			{
				if (TreePoint[i].tree == 0)
				{
					if (direction == 0)
					{
						if ((TreePoint[i].x == (TreePoint[processpoint[1]].x + a)) && (TreePoint[i].y == TreePoint[processpoint[1]].y)) { right = TreePoint[i].z; if (direction == 0) processpoint[2] = i; /*cout << "right=" << right << endl;*/ flaglr = 1; }
					}

					if (direction == 1)
					{
						if ((TreePoint[i].x == (TreePoint[processpoint[1]].x - a)) && (TreePoint[i].y == TreePoint[processpoint[1]].y)) { left = TreePoint[i].z; if (direction == 1) processpoint[0] = i; /*cout << "left=" << left << endl;*/ flaglr = 1; }
					}
					//cout << TreePoint[i].x - (TreePoint[processpoint[1]].x ) << endl;
					//cout << TreePoint[i].x - (TreePoint[processpoint[1]].x ) << endl;
				}
			}
			//if (i >= countR) { direction++; }
			if (((direction == 0) && (TreePoint[processpoint[1]].z <= right)) || ((direction == 1) && (TreePoint[processpoint[1]].z <= left))) //�����������б��������
			{
				if (direction == 0) treeright = processpoint[1];  //��ȡ���Ҳ�countֵ
				if (direction == 1) treelift = processpoint[1];  //��ȡ�����countֵ
																 //left = right = 0;
				TreePoint[i].X = treenumber;
				processpoint[1] = numbermax;
				flag = 1;
				direction++;
				//	cout << direction << endl;
			}
			if (flag == 0)
			{
				if (flaglr == 0)  //��Ѱʧ��
				{
					if (direction == 0) treeright = processpoint[1];  //��ȡ���Ҳ�countֵ
					if (direction == 1) treelift = processpoint[1];  //��ȡ�����countֵ
					direction++;
				}
				if (flaglr == 1)
				{
					if (direction == 0) processpoint[1] = processpoint[2];
					if (direction == 1) processpoint[1] = processpoint[0];
				}
			}
			flag = 0;
			flaglr = 0;
		}
		direction = 2;
		//cout << treeright << "   " << treelift << endl;
		processpoint[2] = processpoint[0] = 0;

		for (i = 0; i < countR; i++)  //y�᷽������
		{
			if (TreePoint[i].tree == 0)
			{
				processpoint[1] = i;
				if ((TreePoint[i].x >= TreePoint[treelift].x) && (TreePoint[i].x <= TreePoint[treeright].x) && (TreePoint[i].y == TreePoint[numbermax].y))
				{
					//TreePoint[i].X = treenumber;
					while ((direction == 2) || (direction == 3))  //���·���
					{
						for (j = 0; j < countR; j++)  //
						{
							//if (TreePoint[j].tree == 0)
							//{
							if (direction == 3)
								if ((TreePoint[j].y == (TreePoint[processpoint[1]].y - a)) && (TreePoint[j].x == TreePoint[processpoint[1]].x)) { down = TreePoint[j].z; if (direction == 3) processpoint[0] = j; /*cout << "down=" << down << endl*/;  break; }
							if (direction == 2)
								if ((TreePoint[j].y == (TreePoint[processpoint[1]].y + a)) && (TreePoint[j].x == TreePoint[processpoint[1]].x)) { up = TreePoint[j].z; if (direction == 2)processpoint[2] = j; /*cout << "up=" << up << endl*/; break; }
							//}
						}
						//	if (j >= countR) direction++;
						if (((direction == 3) && (TreePoint[processpoint[1]].z <= down)) || ((direction == 2) && (TreePoint[processpoint[1]].z <= up))) //�����������б��������
						{
							if (direction == 3) { treedown[k] = processpoint[1]; k++; }//��ȡ���²�countֵ
							if (direction == 2) { treeup[p] = processpoint[1]; p++; } //��ȡ���ϲ�countֵ
																					  //TreePoint[j].X = 1;
							processpoint[1] = i;
							flag = 1;
							direction++;
							//down = up = 0;
							//cout << direction << endl;
						}
						if (flag == 0)
						{
							if (direction == 3) processpoint[1] = processpoint[0];
							if (direction == 2) processpoint[1] = processpoint[2];
						}
						flag = 0;
					}
				}
				direction = 2;
			}
		}
		//cout << "k=" << k << " p=" << p << endl;
		int jj;
		for (jj = 0; jj < k; jj++)   //��ʶ��ĵ�ľ���б��
		{
			for (i = 0; i < countR; i++)
			{
				int kk = treedown[jj];
				int qq = treeup[jj];
				if (TreePoint[i].x == (TreePoint[treelift].x + jj*a))

				{
					if ((TreePoint[i].y >= TreePoint[kk].y) && (TreePoint[i].y <= TreePoint[qq].y))
						/*if (TreePoint[i].y== TreePoint[numbermax].y)*/
					{
						TreePoint[i].X = treenumber;
						TreePoint1[i].X = treenumber;
					}
				}
			}
		}

		diameterlf = TreePoint[treeright].x - TreePoint[treelift].x; //��¼x�᳤��
																	 //TreePoint[numbermax].X = treenumber;
																	 ///*****************��YΪ�����ʶ��*********************/

																	 /*cout << "��ʼ��YΪ��" << endl;*/

		left = 0;
		right = 0;
		up = 0;
		down = 0;
		flag = 0;
		direction = 0;// 0���� 1����  2���� 3����
		processpoint[1] = numbermax;
		k = p = 0;
		while (direction < 2)  //x���ұ߽��
		{
			//if (flag == 1) { processpoint = numbermax; flag = 0; }
			for (i = 0; i < countR; i++)  //��ȡ��ǰ�����¸���Ĵ�С
			{
				if (TreePoint[i].tree == 0)
				{
					if (direction == 1)
						if ((TreePoint[i].y == (TreePoint[processpoint[1]].y - a)) && (TreePoint[i].x == TreePoint[processpoint[1]].x)) { left = TreePoint[i].z; processpoint[0] = i;/* cout << "down=" << left << endl;*/ flaglr = 1; }
					if (direction == 0)
						if ((TreePoint[i].y == (TreePoint[processpoint[1]].y + a)) && (TreePoint[i].x == TreePoint[processpoint[1]].x)) { right = TreePoint[i].z;  processpoint[2] = i; /*cout << "up=" << right << endl;*/ flaglr = 1; }
					//cout << TreePoint[i].x - (TreePoint[processpoint[1]].x ) << endl;
					//cout << TreePoint[i].x - (TreePoint[processpoint[1]].x ) << endl;
				}

			}

			// if (i>= countR) direction++;
			if (((direction == 1) && (TreePoint[processpoint[1]].z <= left)) || ((direction == 0) && (TreePoint[processpoint[1]].z <= right))) //�����������б��������
			{
				if (direction == 0) treeright = processpoint[1];  //��ȡ���Ҳ�countֵ
				if (direction == 1) treelift = processpoint[1];  //��ȡ�����countֵ
																 //left = right = 0;
																 // TreePoint[i].X = treenumber;
				processpoint[1] = numbermax;
				flag = 1;
				direction++;
				// cout << direction << endl;
			}
			if (flag == 0)
			{
				if (flaglr == 0)  //��Ѱʧ��
				{
					if (direction == 0) treeright = processpoint[1];  //��ȡ���Ҳ�countֵ
					if (direction == 1) treelift = processpoint[1];  //��ȡ�����countֵ
					direction++;
				}
				if (flaglr == 1)
				{
					if (direction == 0) processpoint[1] = processpoint[2];
					if (direction == 1) processpoint[1] = processpoint[0];
				}
			}
			flag = 0;
			flaglr = 0;
		}
		direction = 2;
		// cout << treeright << "   " << treelift << endl;
		processpoint[2] = processpoint[0] = 0;

		for (i = 0; i < countR; i++)  //y�᷽������
		{
			if (TreePoint[i].tree == 0)
			{
				processpoint[1] = i;
				if ((TreePoint[i].y >= TreePoint[treelift].y) && (TreePoint[i].y <= TreePoint[treeright].y) && (TreePoint[i].x == TreePoint[numbermax].x))
				{
					//TreePoint[i].X = treenumber;
					while ((direction == 2) || (direction == 3))  //���·���
					{
						for (j = 0; j < countR; j++)  //
						{
							//if (TreePoint[j].tree == 0)
							//{
							if (direction == 3)
								if ((TreePoint[j].x == (TreePoint[processpoint[1]].x - a)) && (TreePoint[j].y == TreePoint[processpoint[1]].y)) { down = TreePoint[j].z; if (direction == 3) processpoint[0] = j;/* cout << "lift=" << down << endl;*/ break; }
							if (direction == 2)
								if ((TreePoint[j].x == (TreePoint[processpoint[1]].x + a)) && (TreePoint[j].y == TreePoint[processpoint[1]].y)) { up = TreePoint[j].z; if (direction == 2)processpoint[2] = j;/* cout << "right=" << up << endl;*/ break; }
							//}
						}
						// if (j >= countR) direction++;
						if (((direction == 3) && (TreePoint[processpoint[1]].z <= down)) || ((direction == 2) && (TreePoint[processpoint[1]].z <= up))) //�����������б��������
						{
							if (direction == 3) { treedown[k] = processpoint[1]; k++; }//��ȡ���²�countֵ
							if (direction == 2) { treeup[p] = processpoint[1]; p++; } //��ȡ���ϲ�countֵ
																					  //TreePoint[j].X = 1;
							processpoint[1] = i;
							flag = 1;
							direction++;
						}
						if (flag == 0)
						{
							if (direction == 3) processpoint[1] = processpoint[0];
							if (direction == 2) processpoint[1] = processpoint[2];
						}
						flag = 0;
					}
				}
				direction = 2;
			}
		}

		int jjj;
		for (jjj = 0; jjj < k; jjj++)   //Y��ʶ��ĵ�ľ���б��
		{
			for (i = 0; i < countR; i++)
			{
				int kk = treedown[jjj];
				int qq = treeup[jjj];
				if (TreePoint[i].y == (TreePoint[treelift].y + jjj*a))

				{
					if ((TreePoint[i].x >= TreePoint[kk].x) && (TreePoint[i].x <= TreePoint[qq].x))
						TreePoint[i].tree = treenumber;
					TreePoint1[i].tree = treenumber;
				}
			}
		}
		diameterud = TreePoint[treeright].y - TreePoint[treelift].y;  //��¼y�᳤��  



		if (diameterud > diameterlf) diametermax = diameterud;
		else diametermax = diameterlf;
		// cout << "diametermax" << diameterlf << endl;
		// cout << "diametermax" << diameterud << endl;
		cout << "diametermax" << diametermax << endl;
		TreePoint[numbermax].tree = treenumber;
		TreePoint[numbermax].maxflag = 1;
		for (i = 0; i < countR; i++)
		{
			if (TreePoint[i].X != 0)
			{
				TreePoint[i].tree = TreePoint[i].X;
				TreePoint1[i].tree = TreePoint1[i].X;
			}
		}

		for (i = 0; i < countR; i++)  //�������X,Y��������ȡ��
		{
			if (TreePoint[i].tree == treenumber)
			{
				if ((pow(TreePoint[numbermax].x - TreePoint[i].x, 2) + pow(TreePoint[numbermax].y - TreePoint[i].y, 2))>(1.2*pow(0.5*diametermax, 2)))
				{
					TreePoint[i].tree = 0;
					TreePoint[i].X = 0;
					TreePoint1[i].tree = 0;
					TreePoint1[i].X = 0;
					// cout <<"ʵ�ʾ���" <<(pow(TreePoint[numbermax].x - TreePoint[i].x, 2) + pow(TreePoint[numbermax].y - TreePoint[i].y, 2)) << endl;
					// cout << "������ֵ"<<(1.2* pow(0.5*diametermax, 2)) << endl;
				}
			}
		}

		double treenumber_sum = 0, treenumber_0 = 0;
		for (i = 0; i < countR; i++)  //�������X,Y��������ȡ��
		{
			if ((pow(TreePoint[numbermax].x - TreePoint[i].x, 2) + pow(TreePoint[numbermax].y - TreePoint[i].y, 2)) <= (1.2*pow(0.5*diametermax, 2)))
			{
				treenumber_sum++;
				if (TreePoint[i].tree != treenumber) treenumber_0++;
			}
		}
		treenumber_sum = treenumber_0 / treenumber_sum;
		for (i = 0; i < countR; i++)  //�������X,Y��������ȡ��
		{
			if (TreePoint[i].tree == treenumber)
				TreePoint[i].Canopy_ratio = treenumber_sum;
			TreePoint[i].crown = 2 * (1.2*pow(0.5*diametermax, 2));
		}

		cout << "�����:" << treenumber << endl;
	}
	return treenumber - 1;
}

void  ScreeningOfTree(int numbertree, int MinNumberOfPoint)  //��ʶ���������ɸѡ �����ܿ�������������������Сֵ������Ϊ5�����������
{
	int i = 0, j = 0, k = 0, p = 0, pointnumber = 0; //��ľ�����
	int maxpoint = 0;
	int countmin = 0; //��¼�����������countֵ
	double mindistance = 9999.9, distance; //��С����͵�ǰ����
	for (i = 1; i < numbertree + 1; i++)
	{
		for (j = 0; j < countR; j++)
		{
			if (TreePoint[j].tree == i)
			{
				pointnumber++;
				if (TreePoint[j].maxflag == 1) maxpoint = j;  //��¼��ǰ����countֵ

			} //��¼��ĸ��������ҵ����ֵ��
		}
		if (pointnumber < MinNumberOfPoint)//������Ƹ���С����С�����������ҵ�������ӽ������ֵ���ϲ�Ϊͬһ����
		{
			for (k = 0; k < countR; k++)
			{
				if ((TreePoint[k].maxflag == 1) && (k != maxpoint)) //Ϊ���㣬�Ҳ��ǵ�ǰ��
				{
					distance = pow((TreePoint[k].x - TreePoint[maxpoint].x), 2) + pow((TreePoint[k].y - TreePoint[maxpoint].y), 2);
					if (mindistance>distance) //�����С����Ϊ��ǰ����,�򽫵��Ϊ����
					{
						mindistance = distance;
						countmin = k;//��¼�����������countֵ
					}
				}
			}
			for (p = 0; p < countR; p++)
			{
				if (TreePoint[p].tree == i)
				{
					TreePoint[p].tree = TreePoint[countmin].tree;
					TreePoint[p].maxflag = 0;
				}
			}
		}
		pointnumber = 0;
	}
}

void ShapeDetection(int numbertree)  //������Բ�����ԣ�����������״ɸѡ
{
	int i, j, pointnumber = 0;
	double meanx = 0, meany = 0;  //x,yƽ��ֵ
	double variancex = 0, variancey = 0;
	double d = 0;
	int n = 0; //��¼�����ܵ���
	double dthreshold;
	dthreshold = 1.55 - 0.5*a;
	for (i = 1; i < numbertree + 1; i++)
	{
		for (j = 0; j < countR; j++)
		{
			if (TreePoint[j].tree == i)
			{
				pointnumber++;
				meanx += 0.1*TreePoint[j].x;
				meany += 0.1*TreePoint[j].y;
			} //��¼��ĸ���

		}
		meanx = 10 * (meanx / pointnumber);
		meany = 10 * (meany / pointnumber);
		//cout << "meanx=" << meanx << endl;
		//cout << "meany=" << meany << endl;
		for (j = 0; j < countR; j++)
		{
			if (TreePoint[j].tree == i)
			{
				variancex += pow((TreePoint[j].x - meanx), 2);
				variancey += pow((TreePoint[j].y - meany), 2);
				n++;
				/*	cout << "meanx=" << variancex << endl;
				cout << "meany=" << variancey << endl;*/
			}
		}

		//cout << "n=" << n << endl;
		//cout << "1+sqrt=" << 1 + sqrt(variancex + variancey) << endl;
		d = (double)n / (1 + sqrt(variancex + variancey)); //�ܶ�����d
														   //cout << "d=" <<d<<endl;
		if ((d< dthreshold) || (n <8))
		{
			for (j = 0; j < countR; j++)
			{
				if (TreePoint[j].tree == i)
				{
					TreePoint[j].tree = 0;
				}
			}
		}
		n = 0;
		pointnumber = 0;
		meanx = 0;
		meany = 0;
		variancex = 0;
		variancey = 0;
	}

	for (i = 0; i < countR; i++)
	{
		if (TreePoint[i].tree != 0)
		{
			PointRedundancy[i].b = 30 + TreePoint[i].tree * 8;
			PointRedundancy[i].r = 50 + TreePoint[i].tree * 10;
			PointRedundancy[i].g = 10 + TreePoint[i].tree * 7;
		}
		else
		{
			PointRedundancy[i].b = 0;
			PointRedundancy[i].r = 0;
			PointRedundancy[i].g = 0;
		}
	}
}

int SortingLabel(int treenumber) //��������ľ��ǩ,������ľ����
{
	int i, j, k;

	bool flag;
	k = 0;
	for (i = 1; i < treenumber + 1; i++)
	{

		for (j = 0; j < countR; j++)
		{
			if (TreePoint[j].tree == i)
			{
				flag = 1;
				break;
			}
		}

		if (flag == 0)
		{
			cout << "kong" << i << endl;
			for (j = 0; j < countR; j++)
			{
				if (TreePoint[j].tree>i)
				{
					TreePoint[j].tree -= 1;

				}
			}
			treenumber -= 1;
			i = i - 1;
		}
		flag = 0;
	}

	return treenumber;
}

int Segmentation(double number, int tree)  //�����ܵ��Ƹ��� ,����� ,���ص�ľ������
{
	long i = 0;

	int pointnumber = 0;
	
			for (i = 0; i < number; i++)
			{
				if (AllPoint[i].tree == tree)
				{
					//AllPoint[i].tree = AllPoint[i].b;
					//AllPoint[i].Canopy_ratio = TreePoint[j].Canopy_ratio;
					//AllPoint[i].crown = TreePoint[j].crown;
					pointnumber++;
				}
			}
	
	return pointnumber;
}
/*�ṹ����*/
void StructuralFeatures(double number, int tree, int pointnumber) //�����ܵ��Ƹ���, ����ţ���ľ������
{
	int i;

	double MaxPointZ = 0;
	double MinPointZ = 10000;
	double Hmed = 0, Hmean = 0, H10, H30, H50, H70, H90, H20, H40, H60, H80, H95;
	double R10 = 0, R30 = 0, R50 = 0, R70 = 0, R90 = 0, R20 = 0, R40 = 0, R60 = 0, R80 = 0, R95 = 0;
	double StandardDeviation = 0;  //��׼��
	double Kurtosis = 0, Skewness = 0; //��ȣ� ƫ��
	double Mean_Intensity = 0, Mean_FirstIntensity = 0;
	int FirstPoint_Num = 0;
	int med_mumber = 0;
	double z_array[50000] = { 0 };
	double temp = 0;
	int tree_count = 0;

	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{
			if (AllPoint[i].z>MaxPointZ) MaxPointZ = AllPoint[i].z;   //�����Zֵ
			if (AllPoint[i].z<MinPointZ) MinPointZ = AllPoint[i].z;
			z_array[tree_count] = AllPoint[i].z;  //��zֵ�浽������
			Hmean += AllPoint[i].z;
			tree_count++;
		}
	}
	med_mumber = tree_count / 2;
	for (i = 0; i<tree_count - 1; i++)
	{
		for (int j = i + 1; j<tree_count; j++)
			if (z_array[j]<z_array[i])
			{
				temp = z_array[i];
				z_array[i] = z_array[j];
				z_array[j] = temp;
			}
	}
	Hmed = z_array[med_mumber];

	Hmean = Hmean / pointnumber;  //��ƽ��ֵ
	H20 = MaxPointZ*0.2;
	H40 = MaxPointZ*0.4;
	H60 = MaxPointZ*0.6;
	H80 = MaxPointZ*0.8;
	H95 = MaxPointZ*0.95;
	H10 = MaxPointZ*0.1;
	H30 = MaxPointZ*0.3;
	H50 = MaxPointZ*0.5;
	H70 = MaxPointZ*0.7;
	H90 = MaxPointZ*0.9;
	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{

			if (AllPoint[i].z <= H10) R10++;
			else if ((AllPoint[i].z > H10) && (AllPoint[i].z <= H20)) R20++;
			else if ((AllPoint[i].z > H20) && (AllPoint[i].z <= H30)) R30++;
			else if ((AllPoint[i].z > H30) && (AllPoint[i].z <= H40)) R40++;
			else if ((AllPoint[i].z > H40) && (AllPoint[i].z <= H50)) R50++;
			else if ((AllPoint[i].z > H50) && (AllPoint[i].z <= H60)) R60++;
			else if ((AllPoint[i].z > H60) && (AllPoint[i].z <= H70)) R70++;
			else if ((AllPoint[i].z > H70) && (AllPoint[i].z <= H80)) R80++;
			else if ((AllPoint[i].z > H80) && (AllPoint[i].z <= H90)) R90++;
			else if ((AllPoint[i].z > H90) && (AllPoint[i].z <= H95)) R95++;
			StandardDeviation += pow(AllPoint[i].z - Hmean, 2);
			Skewness += pow(AllPoint[i].z - Hmean, 3);
			Kurtosis += pow(AllPoint[i].z - Hmean, 4);
			Mean_Intensity += AllPoint[i].intensity;
			if (AllPoint[i].ReturnNumber == 1)
			{
				Mean_FirstIntensity += AllPoint[i].intensity;
				FirstPoint_Num++;
			}
			//Canopy_ratio = AllPoint[i].Canopy_ratio;
		}
	}
	R20 = R20 / pointnumber;
	R40 = R40 / pointnumber;
	R60 = R60 / pointnumber;
	R80 = R80 / pointnumber;
	R95 = R95 / pointnumber;
	R10 = R10 / pointnumber;
	R30 = R30 / pointnumber;
	R50 = R50 / pointnumber;
	R70 = R70 / pointnumber;
	R90 = R90 / pointnumber;

	StandardDeviation = sqrt(StandardDeviation / pointnumber);

	Skewness = (Skewness / pow(StandardDeviation, 3)) / pointnumber;
	Kurtosis = (Kurtosis / pow(StandardDeviation, 4)) / pointnumber;
	Mean_Intensity = Mean_Intensity / pointnumber;  //���е�ǿ�Ⱦ�ֵ
	Mean_FirstIntensity = Mean_FirstIntensity / FirstPoint_Num; //�׻ز�ǿ�Ⱦ�ֵ



	
	MatrixTransform(number, tree, 4, Kurtosis);  //���
	MatrixTransform(number, tree, 5, Skewness);  //ƫ��
	MatrixTransform(number, tree, 6, StandardDeviation);  //��׼��

	MatrixTransform(number, tree, 7, MaxPointZ); //�߶����ֵ
	MatrixTransform(number, tree, 8, Hmed);
	MatrixTransform(number, tree, 9, Hmean);  //�߶Ⱦ�ֵ

	MatrixTransform(number, tree, 10, R10);
	MatrixTransform(number, tree, 11, R20);
	MatrixTransform(number, tree, 12, R30);
	MatrixTransform(number, tree, 13, R40);
	MatrixTransform(number, tree, 14, R50);
	MatrixTransform(number, tree, 15, R60);
	MatrixTransform(number, tree, 16, R70);
	MatrixTransform(number, tree, 17, R80);
	MatrixTransform(number, tree, 18, R90);
	MatrixTransform(number, tree, 19, R95);

}
/*��������*/
void TexturalFeatures(double number, int tree, int pointnumber)//���� ����ţ���ľ������
{
	int i, j = 0, k, q, p, r, s;
	double MaxPointZ = 0, MaxPointX = 0, MaxPointY = 0;
	double MinPointZ = 10000, MinPointX = 10000, MinPointY = 10000;
	int ceilx, ceily, ceilz;
	int maxceil = 0;
	int d1, d2, d3;
	unsigned int pointnumber1 = pointnumber;
	int number1 = number;
	struct TreePoint1
	{
		double x = 0;
		double y = 0;
		double z = 0;
	};
	TreePoint1 *TFpoint = NULL;
	TFpoint = new  TreePoint1[pointnumber1];




	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{
			if (AllPoint[i].x>MaxPointX) MaxPointX = AllPoint[i].x;   //�����Xֵ
			if (AllPoint[i].x<MinPointX) MinPointX = AllPoint[i].x;
			if (AllPoint[i].y>MaxPointY) MaxPointY = AllPoint[i].y;   //�����Yֵ
			if (AllPoint[i].y<MinPointY) MinPointY = AllPoint[i].y;
			if (AllPoint[i].z>MaxPointZ) MaxPointZ = AllPoint[i].z;   //�����Zֵ
			if (AllPoint[i].z<MinPointZ) MinPointZ = AllPoint[i].z;
		}
	}
	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{
			TFpoint[j].x = AllPoint[i].x - MinPointX;
			TFpoint[j].y = AllPoint[i].y - MinPointY;
			TFpoint[j].z = AllPoint[i].z - MinPointZ;
			//cout << "pointx=" << point[j].x;
			j++;
		}
	}
	//cout << "j-pointnumber = " << j - pointnumber << endl;
	MaxPointX = MaxPointX - MinPointX;
	MaxPointY = MaxPointY - MinPointY;
	MaxPointZ = MaxPointZ - MinPointZ;

	ceilx = MaxPointX / 0.5;
	ceily = MaxPointY / 0.5;
	ceilz = MaxPointZ / 0.5;


	int *** voxel;  //�������ؾ��� ��¼�����
	voxel = new int**[ceilx + 1];
	for (i = 0; i <= ceilx; i++)
		voxel[i] = new int*[ceily + 1];
	for (i = 0; i <= ceilx; i++)
		for (j = 0; j <= ceily; j++)
			voxel[i][j] = new int[ceilz + 1];

	for (i = 0; i <= ceilx; i++)
		for (j = 0; j <= ceily; j++)
			for (k = 0; k <= ceilz; k++)
				voxel[i][j][k] = 0;
	int ceilx1;
	int ceily1;
	int ceilz1;
	for (i = 0; i < pointnumber; i++)    //ȷ��ÿ�������е�������
	{
		ceilx1 = TFpoint[i].x / 0.5;
		ceily1 = TFpoint[i].y / 0.5;
		ceilz1 = TFpoint[i].z / 0.5;
		voxel[ceilx1][ceily1][ceilz1]++;

	}




	for (i = 0; i <= ceilx; i++)  //ȷ���Ҷȼ�
		for (j = 0; j <= ceily; j++)
			for (k = 0; k <= ceilz; k++)
				if (voxel[i][j][k]>maxceil)  maxceil = voxel[i][j][k];

	int *** GLCM;  //�����Ҷȹ�������
	GLCM = new int**[maxceil];
	for (i = 0; i < maxceil; i++)
		GLCM[i] = new int*[maxceil];
	for (i = 0; i < maxceil; i++)
		for (j = 0; j < maxceil; j++)
			GLCM[i][j] = new int[13];

	for (i = 0; i < maxceil; i++)
		for (j = 0; j < maxceil; j++)
			for (k = 0; k <13; k++)
				GLCM[i][j][k] = 0;


	for (int w = 0; w < 13; w++)
	{
		switch (w)
		{
		case 0: d1 = 1; d2 = 0; d3 = 0; break;
		case 1: d1 = 0; d2 = 1; d3 = 0; break;
		case 2: d1 = 1; d2 = 1; d3 = 0; break;
		case 3: d1 = 1; d2 = -1; d3 = 0; break;
		case 4: d1 = 0; d2 = 0; d3 = 1; break;
		case 5: d1 = 0; d2 = 1; d3 = 1; break;
		case 6: d1 = 0; d2 = 1; d3 = -1; break;
		case 7: d1 = 1; d2 = 0; d3 = -1; break;
		case 8: d1 = 1; d2 = 1; d3 = 1; break;
		case 9: d1 = 1; d2 = -1; d3 = -1; break;
		case 10: d1 = 1; d2 = 1; d3 = 1; break;
		case 11: d1 = 1; d2 = 1; d3 = -1; break;
		case 12: d1 = 1; d2 = -1; d3 = 1; break;
		}
		for (i = 1; i <= maxceil; i++)
			for (j = 1; j <= maxceil; j++)
			{
				//if ((i != 0) || (j != 0))
				//{
				for (k = 0; k <= ceilx; k++)
					for (q = 0; q <= ceily; q++)
						for (p = 0; p <= ceilz; p++)
						{
							if (((k + d1 >= 0) && (k + d1 <= ceilx)) && ((q + d2 >= 0) && (q + d2 <= ceily)) && ((p + d3 >= 0) && (p + d3 <= ceilz)))
							{
								if ((voxel[k][q][p] == i) && (voxel[k + d1][q + d2][p + d3] == j))
								{
									GLCM[i - 1][j - 1][w]++;
								}
							}
						}
				//}

			}
	}
	double f1 = 0;//�Ƕ��׾�
	double f2 = 0;
	double u1 = 0, u2 = 0, q1 = 0, q2 = 0, f31 = 0, f3 = 0;
	double f4 = 0;
	double f5 = 0, f6 = 0, f7 = 0, f8 = 0, f9 = 0, f10 = 0, f11 = 0, f12 = 0, f13 = 0, f14 = 0;
	double Pxy = 0, Pxy0 = 0, Px = 0, Py = 0;
	double DIFave = 0, HX = 0, HY = 0, HXY1 = 0, HXY2 = 0;
	//cout << "�Ҷȹ�������" << endl;
	for (int w = 0; w < 13; w++)
	{
		//cout << "����" << w << endl;

		for (i = 1; i <= maxceil; i++)
			for (j = 1; j <= maxceil; j++)
			{
				u1 += i*GLCM[i - 1][j - 1][w];
				u2 += j*GLCM[i - 1][j - 1][w];
			}
		u1 = u1 / maxceil;
		u2 = u2 / maxceil;
		for (i = 1; i <= maxceil; i++)
			for (j = 1; j <= maxceil; j++)
			{
				q1 += GLCM[i - 1][j - 1][w] * pow((i - u1), 2);
				q2 += GLCM[i - 1][j - 1][w] * pow((j - u2), 2);
			}
		q1 = sqrt(q1);
		q2 = sqrt(q2);

		for (i = 1; i <= maxceil; i++)
		{
			for (j = 1; j <= maxceil; j++)
			{
				//cout << GLCM[i-1][j-1][w] << " ";
				f1 += pow(GLCM[i - 1][j - 1][w], 2);
				for (int n = 0; n < maxceil - 1; n++)
				{
					if (pow(i - j, 2) == pow(n, 2))
						f2 += pow(i - j, 2)*GLCM[i - 1][j - 1][w];
				}
				f31 += i*j*GLCM[i - 1][j - 1][w];
				f4 += pow((i - u1), 2)*GLCM[i - 1][j - 1][w];
				f5 += (1 / (1 + pow(i + j, 2)))*GLCM[i - 1][j - 1][w];
				for (int k = 2; k <= 2 * maxceil; k += 2)
				{
					if (i + j == k) f6 += k*GLCM[i - 1][j - 1][w];
				}
				f9 += GLCM[i - 1][j - 1][w] * log(GLCM[i - 1][j - 1][w] + 1);
			}
		}
		for (int k = 2; k <= 2 * maxceil; k += 2)
		{
			Pxy = 0;
			for (i = 1; i <= maxceil; i++)
			{
				for (j = 1; j <= maxceil; j++)
				{
					if (i + j == k) Pxy += GLCM[i - 1][j - 1][w];
				}
			}
			f10 += Pxy*log(Pxy + 1);

		}

		f31 = (f31 - u1*u2*maxceil*maxceil) / (q1*q2);
		f3 += f31;
		f31 = 0;
	}
	f1 = f1 / 13.0;
	f2 = f2 / 13.0;
	f3 = f3 / 13.0;
	f4 = f4 / 13.0;
	f5 = f5 / 13.0;
	f6 = f6 / 13.0;
	f9 = f9 / 13.0;
	f10 = -f10 / 13.0;

	for (int w = 0; w < 13; w++)
	{
		for (int k = 2; k <= 2 * maxceil; k += 2)
		{
			Pxy = 0;
			for (i = 1; i <= maxceil; i++)
			{
				for (j = 1; j <= maxceil; j++)
				{
					if (i + j == k) Pxy += GLCM[i - 1][j - 1][w];
				}
			}
			f7 += pow(k - f6, 2)*Pxy;
		}

		for (int k = 0; k <= maxceil - 1; k++)
		{
			Pxy0 = 0;
			for (i = 1; i <= maxceil; i++)
			{
				for (j = 1; j <= maxceil; j++)
				{
					if (pow(i - j, 2) == pow(k, 2)) Pxy0 += GLCM[i - 1][j - 1][w];
				}
			}
			DIFave += k*Pxy0;
		}

		for (int k = 0; k <= maxceil - 1; k++)
		{
			Pxy0 = 0;
			for (i = 1; i <= maxceil; i++)
			{
				for (j = 1; j <= maxceil; j++)
				{
					if (pow(i - j, 2) == pow(k, 2)) Pxy0 += GLCM[i - 1][j - 1][w];
				}
			}
			f8 += pow(k - DIFave, 2)*Pxy0;
			f11 += Pxy0*log(Pxy0 + 1);
		}
		for (i = 1; i <= maxceil; i++)
		{
			Px = 0;
			for (j = 1; j <= maxceil; j++)
			{
				Px += GLCM[i - 1][j - 1][w];
			}

			for (r = 1; r <= maxceil; r++)
			{
				Py = 0;
				for (s = 1; s <= maxceil; s++)
				{
					Py += GLCM[s - 1][r - 1][w];
				}
				HY += Py*log(Py + 1) / maxceil;
				HXY1 += GLCM[i - 1][r - 1][w] * log(Px*Py + 1) / maxceil;
				HXY2 += Px*Py*log(Px*Py + 1) / maxceil;
				
			}
			HX += Px*log(Px + 1);
		}
		if (HY > HX) HX = HY;
		f12 += (f9 - HXY1) / HX;
		HX = 0;
		HY = 0;
		Px = 0;
		Py = 0;
		f13 += pow(1 - exp(-2.0 * (-HXY2 - f9)), 0.5);
		

		for (k = 1; k <= maxceil; k++)
		{
			for (i = 1; i <= maxceil; i++)
			{
				for (j = 1; j <= maxceil; j++)
				{
					s = GLCM[i - 1][k - 1][w] * GLCM[j - 1][k - 1][w];

					Px = 0;
					Py = 0;
					for (r = 1; r <= maxceil; r++)
					{
						Px += GLCM[i - 1][r - 1][w];
						Py += GLCM[r - 1][k - 1][w];
					}
					if (Px*Py != 0) f14 += s / (Px*Py) / maxceil;
				}
			}
		}
	}

	f7 = f7 / 13.0;
	f8 = f8 / 13.0;
	f11 = -f11 / 13.0;
	f12 = f12 / 13.0;
	f13 = f13 / 13.0;
	f14 = f14 / 13.0;
	//cout <<  "�Ƕ��׾�f1 = " << f1 << endl;
	//cout <<  "  �Աȶ�f2 = " << f2 << endl;
	//cout << "   �����f3 = " << f3 << endl;
	//cout << "ƽ���ͷ���f4 = " << f4 << endl;
	cout << "   ����f5 = " << HXY2 << endl;
	
	MatrixTransform(number, tree, 20, f1);
	MatrixTransform(number, tree, 21, f2);
	MatrixTransform(number, tree, 22, f3);
	MatrixTransform(number, tree, 23, f4);
	MatrixTransform(number, tree, 24, f5);
	MatrixTransform(number, tree, 25, f6);
	MatrixTransform(number, tree, 26, f7);
	MatrixTransform(number, tree, 27, f8);
	MatrixTransform(number, tree, 28, f9);
	MatrixTransform(number, tree, 29, f10);
	MatrixTransform(number, tree, 30, f11);
	MatrixTransform(number, tree, 31, f12);
	MatrixTransform(number, tree, 32, f13);
	MatrixTransform(number, tree, 33, f14);

	//MatrixTransform(number, tree, 11, f11);
}
/*��������*/
void CanopyFeatures(double number, int tree, int pointnumber)
{
	int i, j = 0, k, p, q;
	double MaxPointZ = 0, MaxPointX = 0, MaxPointY = 0;
	double MinPointZ = 10000, MinPointX = 10000, MinPointY = 10000;
	double w_x = 0, w_y = 0, w_z = 0;
	int ceilx, ceily, ceilz;

	unsigned int pointnumber1 = pointnumber;
	int number1 = number;

	double CF1 = 0, CF2 = 0, CF3 = 0, CF4 = 0, CF5 = 0, CF6 = 0, CF7 = 0, CF8 = 0, CF9 = 0, CF10 = 0, CF11 = 0, CF12 = 0;
	double crown = 0;
	int ceilx1;
	int ceily1;
	int ceilz1;


	struct TreePoint1
	{
		double x = 0;
		double y = 0;
		double z = 0;
	};
	TreePoint1 *CFpoint = NULL;
	CFpoint = new  TreePoint1[pointnumber1];




	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{
			if (AllPoint[i].x>MaxPointX) MaxPointX = AllPoint[i].x;   //�����Xֵ
			if (AllPoint[i].x<MinPointX) MinPointX = AllPoint[i].x;
			if (AllPoint[i].y>MaxPointY) MaxPointY = AllPoint[i].y;   //�����Yֵ
			if (AllPoint[i].y<MinPointY) MinPointY = AllPoint[i].y;
			if (AllPoint[i].z>MaxPointZ) MaxPointZ = AllPoint[i].z;   //�����Zֵ
			if (AllPoint[i].z<MinPointZ) MinPointZ = AllPoint[i].z;
			crown = AllPoint[i].crown;  //����ڷ�
		}
	}
	for (i = 0; i < number; i++)
	{
		if (AllPoint[i].tree == tree)
		{
			CFpoint[j].x = AllPoint[i].x - MinPointX;
			CFpoint[j].y = AllPoint[i].y - MinPointY;
			CFpoint[j].z = AllPoint[i].z - MinPointZ;
			//cout << "pointx=" << point[j].x;
			j++;
		}
	}
	//cout << "j-pointnumber = " << j - pointnumber << endl;
	MaxPointX = MaxPointX - MinPointX;
	MaxPointY = MaxPointY - MinPointY;
	MaxPointZ = MaxPointZ - MinPointZ;

	crown = sqrt(pow(MaxPointX,2)+ pow(MaxPointY, 2));
	double rx = MaxPointX / 2.0;
	double ry = MaxPointY / 2.0;
	w_x = MaxPointX / 19.999;
	w_y = MaxPointY / 19.999;
	w_z = MaxPointZ / 19.999;

	ceilx = 19;//20*20*20  ������ŵ�19
	ceily = 19;
	ceilz = 19;


	int **** CFvoxel;  //�������ؾ��� ��¼�����
	CFvoxel = new int***[ceilx + 1];
	for (i = 0; i <= ceilx; i++)
		CFvoxel[i] = new int**[ceily + 1];
	for (i = 0; i <= ceilx; i++)
		for (j = 0; j <= ceily; j++)
			CFvoxel[i][j] = new int*[ceilz + 1];
	for (i = 0; i <= ceilx; i++)
		for (j = 0; j <= ceily; j++)
			for (k = 0; k <= ceilz; k++)
				CFvoxel[i][j][k] = new int[2];

	for (i = 0; i <= ceilx; i++)
		for (j = 0; j <= ceily; j++)
			for (k = 0; k <= ceilz; k++)
				for (p = 0; p < 2; p++)
					CFvoxel[i][j][k][p] = 0;
	/*
	voxel[i][j][k][p]
	i x����
	j y����
	k z����
	p 0���ص��Ƹ��� 2 ��t���ֶ�

	*/



	for (i = 0; i < pointnumber; i++)    //ȷ��ÿ�������е�������
	{
		ceilx1 = CFpoint[i].x / w_x;
		ceily1 = CFpoint[i].y / w_y;
		ceilz1 = CFpoint[i].z / w_z;

		CFvoxel[ceilx1][ceily1][ceilz1][0] = CFvoxel[ceilx1][ceily1][ceilz1][0] + 1;
	}

	for (i = 0; i < 20; i++)
		for (j = 0; j < 20; j++)
		{
			if ((0 <= i < 10) && (0 <= j<10) && (i >= j)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 1;
			if ((0 <= i < 10) && (0 <= j<10) && (i<j)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 8;
			if ((10 <= i < 20) && (0 <= j<10) && (j<19 - i)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 2;
			if ((10 <= i < 20) && (0 <= j<10) && (j >= 19 - i)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 3;
			if ((10 <= i < 20) && (10 <= j<20) && (i >= j)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 4;
			if ((10 <= i < 20) && (10 <= j<20) && (i<j)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 5;
			if ((0 <= i < 10) && (10 <= j<20) && (j >= 19 - i)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 6;
			if ((0 <= i < 10) && (10 <= j<20) && (j<19 - i)) for (k = 0; k < 20; k++) CFvoxel[i][j][k][1] = 7;
		}
	double vox_meanz = 0, vox_meanr = 0;
	int vox_number = 0;  //��¼ÿ�����ǿ����ظ���
	double meanz = 0, meanr = 0;
	double number_sumvox = 0;  //���з�0��������
	double number_sumvox31 = 0; //����1/3���Ϸǿ���������
	double number_maxvox = 0; //��¼ÿ�����������ĵ�����
	double meanz_t[8] = { 0 }, meanr_t[8] = { 0 };  //��¼ÿ��
	double voxnum_t[8] = { 0 };  //��¼ÿ�����ǿ����ظ���
	double sub_vox = 0;  //�������ص����������
	double sub_t = 0;  //�������ǿ����ظ������
	for (p = 1; p <= 8; p++)
	{
		for (i = 0; i < 20; i++)
		{
			for (j = 0; j < 20; j++)
			{
				for (k = 0; k < 20; k++)
				{
					if ((CFvoxel[i][j][k][0] != 0) && (CFvoxel[i][j][k][1] == p))
					{
						for (q = 0; q < pointnumber; q++)    //ȷ��ÿ�������е�������
						{
							if ((i == int(CFpoint[q].x / w_x)) && (j == int(CFpoint[q].y / w_y)) && (k == int(CFpoint[q].z / w_z)))
							{
								vox_meanz += CFpoint[q].z;
								vox_meanr += sqrt(pow(CFpoint[q].x - rx, 2) + pow(CFpoint[q].y - ry, 2));
							}
						}
						vox_meanz = vox_meanz / CFvoxel[i][j][k][0];  //ÿ�����ؾ�ֵ		
						vox_meanr = vox_meanr / CFvoxel[i][j][k][0];  //ÿ�����ص�Ч���İ뾶
						meanz += vox_meanz;
						meanr += vox_meanr;
						vox_meanz = 0;
						vox_meanr = 0;
						vox_number++;   //ÿ�����ǿ���������
						number_sumvox++;  //ȫ���ǿ����ظ���
						if (k > 6) number_sumvox31++;   //����1/3���ظ���

						if (CFvoxel[i][j][k][0]>number_maxvox) number_maxvox = CFvoxel[i][j][k][0];
					}
				}
			}
		}

		if (vox_number != 0) { meanz = meanz / vox_number; meanr = meanr / vox_number; }
		else { meanz = 0; meanr = 0; }
		voxnum_t[p - 1] = vox_number;
		meanz_t[p - 1] = meanz / crown;
		meanr_t[p - 1] = meanr / crown * 2;
		CF1 += meanz;
		CF2 += meanr;
		meanz = 0;
		meanr = 0;
		vox_number = 0;
	}

	CF1 += CF1 / (8 * crown);
	CF2 += CF2 / (4 * crown);
	CF3 = number_sumvox31 / number_sumvox;
	CF4 = number_maxvox / (w_x*w_y*w_z);

	for (i = 1; i < 8; i++)
	{
		CF5 += pow(meanz_t[i] - CF1, 2);
	}
	CF5 = sqrt(CF5 / 7);
	for (i = 1; i < 8; i++)
	{
		CF6 += pow(meanr_t[i] - CF2, 2);
	}
	CF6 = sqrt(CF6 / 7);

	double sum_cf8 = 0;
	for (i = 0; i < 8; i++)
	{
		sum_cf8 += voxnum_t[j];
	}
	sum_cf8 = sum_cf8 / 8;
	for (i = 0; i < 8; i++)
	{
		CF8 += pow(voxnum_t[i] - sum_cf8, 2);
	}
	CF8 = sqrt(CF8 / 7);

	for (i = 1; i < 8; i++)
	{
		sub_t += abs(voxnum_t[i] - voxnum_t[i - 1]);
	}
	sub_t += abs(voxnum_t[0] - voxnum_t[7]);
	for (i = 0; i < 20; i++)
	{
		for (j = 0; j < 20; j++)
		{
			for (k = 1; k < 20; k++)
			{
				if ((CFvoxel[i][j][k][0] != 0) && (CFvoxel[i][j][k - 1][0] != 0))
					sub_vox += abs(CFvoxel[i][j][k][0] - CFvoxel[i][j][k - 1][0]);
			}
		}
	}
	CF9 = sub_t / sub_vox;

	double crown_area = 0;
	double crown_volume = 0;
	for (i = 0; i < 20; i++)
	{
		for (j = 0; j < 20; j++)
		{
			for (k = 0; k < 20; k++)
			{
				if (CFvoxel[i][j][k][0] != 0) crown_area++; break;
			}
		}
	}
	for (i = 0; i < 20; i++)
	{
		for (j = 0; j < 20; j++)
		{
			for (k = 0; k < 20; k++)
			{
				if (CFvoxel[i][j][k][0] != 0) crown_volume++;
			}
		}
	}
	CF10 = crown_area / crown_volume;

	CF12 = CF5 / CF6;
	cout << "cf1 = " << CF1 << endl;
	cout << "cf2 = " << CF2 << endl;
	cout << "cf3 = " << CF3 << endl;
	cout << "cf4 = " << CF4 << endl;
	cout << "cf5 = " << CF5 << endl;
	cout << "cf6 = " << CF6 << endl;
	cout << "cf8 = " << CF8 << endl;
	cout << "cf9 = " << CF9 << endl;
	cout << "cf10 = " << CF10 << endl;
	cout << "cf12 = " << CF12 << endl;

	MatrixTransform(number, tree, 34, CF1);
	MatrixTransform(number, tree, 35, CF2);
	MatrixTransform(number, tree, 36, CF3);
	MatrixTransform(number, tree, 37, CF4);
	MatrixTransform(number, tree, 38, CF5);
	MatrixTransform(number, tree, 39, CF6);
	MatrixTransform(number, tree, 40, CF8);
	MatrixTransform(number, tree, 41, CF9);
	MatrixTransform(number, tree, 42, CF10);
	MatrixTransform(number, tree, 43, CF12);
}


void writePointCloudFromLas(const char* strInputLasName, const char* strOutPutPointCloudName, const char* strOutputFeatures)
{

	int NumberOfTree = 0;
	//��las�ļ�
	std::ifstream ifs;
	ifs.open(strInputLasName, std::ios::in | std::ios::binary);
	if (!ifs.is_open())
	{
		std::cout << "�޷���.las" << std::endl;
		return;
	}
	liblas::ReaderFactory readerFactory;
	liblas::Reader reader = readerFactory.CreateWithStream(ifs);
	//д����
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudOutput(new pcl::PointCloud<pcl::PointXYZRGB>());
	cloudOutput->clear();

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr TreePointcloud(new pcl::PointCloud<pcl::PointXYZRGB>());
	TreePointcloud->clear();

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr AllPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>());
	AllPointcloud->clear();

	ofstream fout("E:/LIDAR/project/test40/data2.txt");     //����һ��data.txt���ļ�
																   //double number = reader.GetHeader().GetPointRecordsCount();
																   //unsigned int  number1 = number;
																   //AllPoint = new PointCloud[number1];
																   //for (i = 0; i < number; i++)
																   //AllPoint.
	while (reader.ReadNextPoint())
	{
		double x = reader.GetPoint().GetX();
		double y = reader.GetPoint().GetY();
		double z = reader.GetPoint().GetZ();
		int intensity = reader.GetPoint().GetIntensity();
		uint16_t red = reader.GetPoint().GetColor()[0];
		uint16_t green = reader.GetPoint().GetColor()[1];
		uint16_t blue = reader.GetPoint().GetColor()[2];
		uint16_t  ReturnNumber = reader.GetPoint().GetReturnNumber();
		uint16_t NumberOfReturn = reader.GetPoint().GetNumberOfReturns();
		//if(z>=0.2){
		AllPoint[Allcount].x = x;
		AllPoint[Allcount].y = y;
		AllPoint[Allcount].z = z;
		AllPoint[Allcount].r = red;
		AllPoint[Allcount].g = green;
		AllPoint[Allcount].b = blue;
		AllPoint[Allcount].intensity = intensity;
		AllPoint[Allcount].ReturnNumber = ReturnNumber;
		AllPoint[Allcount].NumberOfReturn = NumberOfReturn;
		Allcount++;

		if (ReturnNumber == 1)
		{
			Point[count].x = x;
			Point[count].y = y;
			Point[count].z = z;
			Point[count].r = red;
			Point[count].g = green;
			Point[count].b = blue;
			Point[count].intensity = intensity;
			Point[count].ReturnNumber = ReturnNumber;
			Point[count].NumberOfReturn = NumberOfReturn;
			fout << count << "  x=" << x << "y = " << y << " z =  " << z << "intensity = " << intensity
				<< " red =" << red << " green =" << green << " blue =" << blue
				<< "  ReturnNumber = " << ReturnNumber << "   NumberOfReturn = " << NumberOfReturn << endl;
			count++;
		}

		//}
	}
	double number = Allcount;
	fout.close();                  //�ر��ļ�
	cout << "point =" << count << endl;
	maxx = max(1, count);  //��������С����
	maxy = max(2, count);
	maxz = max(3, count);
	minx = min(1, count);
	miny = min(2, count);
	minz = min(3, count);
	cout << maxx << " " << minx << endl;
	cout << maxy << " " << miny << endl;
	cout << maxz << " " << minz << endl;

	/*�������*/
	for (int i = 0; i < count; i++)
	{
		Point[i].x -= minx;
		Point[i].y -= miny;
		Point[i].z -= minz;
	}

	for (int i = 0; i < number; i++)
	{
		AllPoint[i].x -= minx;
		AllPoint[i].y -= miny;
		AllPoint[i].z -= minz;
		if (AllPoint[i].x < 0) AllPoint[i].x = 0;
		if (AllPoint[i].y < 0) AllPoint[i].y = 0;
		if (AllPoint[i].z < 0) AllPoint[i].z = 0;
	}

	maxx = max(1, count);  //��������С����
	maxy = max(2, count);
	maxz = max(3, count);
	minx = min(1, count);
	miny = min(2, count);
	minz = min(3, count);

	cout << maxx << " " << minx << endl;
	cout << maxy << " " << miny << endl;
	cout << maxz << " " << minz << endl;
	MultipleX = 39.0 / maxx;
	MultipleY = 39.0 / maxy;

	effective = count;
	for (int i = 0; i < count; i++)
	{
		if (Point[i].z <= 5) effective--;
	}


	//a = sqrt((float)effective /((float)count*((float)count / (float)number)));/*****************������**********************/
	a = (float)sqrt(effective / (count*(count / number)))*0.8;
	cout << "a = " << a << endl;
	a = 1.14;


	int treenumber = 0;

	for (int i = 0; i < number; i++)
	{
		if (treenumber < AllPoint[i].b)  treenumber = AllPoint[i].b;
		AllPoint[i].tree = AllPoint[i].b;
	}

	for (int i = 1; i < treenumber; i++)   //�������ڵ�ֱ��
	{
		double min_x = 9999999, min_y = 9999999,max_x = 0,max_y = 0,r = 0;
		for (int j = 0; j < number; j++)
		{
			if (AllPoint[j].tree == i)
			{
				if (AllPoint[j].x <= min_x)  min_x = AllPoint[j].x;
				if (AllPoint[j].x >= max_x)  max_x = AllPoint[j].x;
				if (AllPoint[j].y <= min_y)  min_y = AllPoint[j].y;
				if (AllPoint[j].y >= max_y)  max_y = AllPoint[j].y;
			}
		}

		r = ((max_y - min_y) + (max_x - min_x)) / 1.8;
		cout << "r = "<<r << endl;
		for (int j = 0; j < number; j++)
		{
			if (AllPoint[j].tree == i)
			{
				AllPoint[j].crown = r;
			}
		}
	}


	H_Matrix(number);
	int singletree = 0;
	while (singletree<treenumber)
	{
		singletree++;
		int Singlepointnumber = Segmentation(number, singletree);  //���ص�ǰ�����Ƹ���
		StructuralFeatures(number, singletree, Singlepointnumber);
		cout << "������������" << Singlepointnumber << endl;
		TexturalFeatures(number, singletree, Singlepointnumber); //����Ҷȹ�������
		CanopyFeatures(number, singletree, Singlepointnumber);
		for (long i = 0; i < number; i++)
		{
			if (AllPoint[i].tree == singletree)
			{
				pcl::PointXYZRGB thePt;  //int rgba = 255<<24 | ((int)r) << 16 | ((int)g) << 8 | ((int)b);
				thePt.x = AllPoint[i].x - 100; thePt.y = AllPoint[i].y; thePt.z = AllPoint[i].z;
				uint32_t rgb = (static_cast<uint32_t>(AllPoint[i].r) << 16 | static_cast<uint32_t>(AllPoint[i].g) << 8 | static_cast<uint32_t>(AllPoint[i].b));
				thePt.rgb = *reinterpret_cast<float*>(&rgb);
				TreePointcloud->push_back(thePt);
				AllPointcloud->push_back(thePt);
			}
		}

		for (long i = 0; i < number; i++)
		{
			pcl::PointXYZRGB thePt;  //int rgba = 255<<24 | ((int)r) << 16 | ((int)g) << 8 | ((int)b);
			thePt.x = AllPoint[i].x; thePt.y = AllPoint[i].y; thePt.z = AllPoint[i].z;
			if (AllPoint[i].tree == singletree) {
				uint32_t rgb = (static_cast<uint32_t>(0) << 16 | static_cast<uint32_t>(0) << 8 | static_cast<uint32_t>(255));
				thePt.rgb = *reinterpret_cast<float*>(&rgb);
			}
			AllPointcloud->push_back(thePt);

		}

	}

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer1(new pcl::visualization::PCLVisualizer("������ʾ"));
	viewer1->setBackgroundColor(1, 1, 1);
	pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb1(AllPointcloud);
	viewer1->addPointCloud<pcl::PointXYZRGB>(AllPointcloud, rgb1, "sample cloud");
	viewer1->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "sample cloud");
	while (!viewer1->wasStopped()) {
		viewer1->spinOnce();
	}
	AllPointcloud->clear();

	cout << countR << endl;
	cout << "����������" << number << endl;
	//cout <<  Allcount << endl;
	cout << "�״λز�����" << count << endl;
	cout << "�״λز�����5�ĵ�����" << effective << endl;
	cout << "�߳���" << a << endl;

	ofstream matout(strOutputFeatures);     //����һ��data.txt���ļ�
	for (int j = 0; j < 40; j++)
	{
		for (int i = 0; i < 40; i++)
		{
			for (int k = 0; k < 44; k++)
			{
				if((lidarmat[j][i][k]<DBL_MAX)&&(lidarmat[j][i][k]>-DBL_MAX))  matout << " " << lidarmat[j][i][k];
				else  matout << " " << 0;
			}
			matout << endl;
		}
		matout << endl;
	}
	matout.close();                  //�ر��ļ�

}



int main(int argc, char *argv[])
{
	//char strInputLasName[] = "C:/Users/25297/Desktop/LiDAR/ABBY1/ABBY_001.las";
	char strOutPutPointCloudName[] = "E:/LIDAR/project/FeatureOfTree/las2.pcd";
	char strInputLasName[255];
	char strOutputFeatures[255];
	strcpy(strInputLasName, argv[1]);
	strcpy(strOutputFeatures, argv[1]);
	strcat(strOutputFeatures, ".txt");
	writePointCloudFromLas(strInputLasName, strOutPutPointCloudName, strOutputFeatures);

	system("pause");
	return 0;
}

