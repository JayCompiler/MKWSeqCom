package seq.comp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import Jama.Matrix;

/*
 * 算法步骤:
 * 1)将原始数据按列组成n行m列矩阵X
 * 2)特征中心化。即每一维的数据都减去该维的均值，使每一维的均值都为0
 * 3)求出协方差矩阵
 * 4)求出协方差矩阵的特征值及对应的特征向量
 * 5)将特征向量按对应的特征值大小从上往下按行排列成矩阵，取前k行组成矩阵p
 * 6)Y=PX 即为降维到k维后的数据
 */
public class PCA {

	private static double threshold = 0.95;// false 特征值阈值
	private static int colsize = 100;

	public double getThreshold() {
		return threshold;
	}

	public void setThreshold(double threshold1) {
		threshold = threshold1;
	}

	public static void main(String[] args) {
		double [] a=new double[3];
		for (int i = 0; i < a.length; i++) {
			System.out.println(a[i]);
		}
	}

	/**
	 * 合并数据
	 * 
	 * @param data1
	 * @param data2
	 * @return
	 */
	public static double[][] mergeData(double[][] data1, double[][] data2) {

		if (data1 == null && data2 != null) {
			return data2;
		} else if (data2 == null && data1 != null) {
			return data1;
		} else if (data1.length == 0 && data2.length != 0) {
			return data2;
		} else if (data1.length != 0 && data2.length == 0) {
			return data1;
		} else if (data1.length != 0 && data2.length != 0) {
			double[][] tdata1 = transData(data1);
			double[][] tdata2 = transData(data2);
			int length1 = tdata1.length;
			int length2 = tdata2.length;
			int row = length1 + length2;
			// System.out.println(tdata1[0].length);
			int col = tdata1[0].length + tdata2[0].length;
			double[][] newdata1 = new double[row][col];
			int count = 0;
			for (int i = 0; i < length1; i++) {
				newdata1[i] = tdata1[i];
				count++;
			}
			for (int i = 0; i < length2; i++) {
				newdata1[count] = tdata2[i];
				count++;
			}
			double[][] newdata = transData(newdata1);
			return newdata;
		} else {
			return null;
		}

	}

	/**
	 * 转置二维数组
	 * 
	 * @param data
	 * @return
	 */
	public static double[][] transData(double[][] data) {
		if (data == null)
			return null;
		int row = data.length;
		int col = data[0].length;
		double[][] newdata = new double[col][row];
		for (int i = 0; i < col; i++)
			for (int j = 0; j < row; j++)
				newdata[i][j] = data[j][i];
		return newdata;
	}

	/**
	 * 
	 * 使每个样本的均值为0
	 * 
	 * @param primary
	 *            原始二维数组矩阵
	 * @return averageArray 中心化后的矩阵
	 */
	public static double[][] changeAverageToZero(double[][] primary) {
		// 行
		int n = primary.length;
		// 列
		int m = primary[0].length;
		double[] sum = new double[m];
		double[] average = new double[m];
		double[][] averageArray = new double[n][m];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				sum[i] += primary[j][i];
			}
			average[i] = sum[i] / n;
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				averageArray[j][i] = primary[j][i] - average[i];
			}
		}
		return averageArray;
	}

	/**
	 * 
	 * 计算协方差矩阵
	 * 
	 * @param matrix
	 *            中心化后的矩阵
	 * @return result 协方差矩阵
	 */
	public static double[][] getVarianceMatrix(double[][] matrix) {
		int n = matrix.length;// 行数
		int m = matrix[0].length;// 列数
		double[][] result = new double[m][m];// 协方差矩阵
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				double temp = 0;
				for (int k = 0; k < n; k++) {
					temp += matrix[k][i] * matrix[k][j];
				}
				result[i][j] = temp / (n - 1);
			}
		}
		return result;
	}

	/**
	 * 求特征值矩阵
	 * 
	 * @param matrix
	 *            协方差矩阵
	 * @return result 向量的特征值二维数组矩阵
	 */
	public static double[][] getEigenvalueMatrix(double[][] matrix) {
		Matrix A = new Matrix(matrix);
		// 由特征值组成的对角矩阵,eig()获取特征值
		// A.eig().getD().print(10, 6);//打印格式，空格，小数点后几位
		double[][] result = A.eig().getD().getArray();
		return result;
	}

	/**
	 * 标准化矩阵（特征向量矩阵）
	 * 
	 * @param matrix
	 *            协方差矩阵
	 * @return result 标准化后的二维数组矩阵
	 */
	public static double[][] getEigenVectorMatrix(double[][] matrix) {
		Matrix A = new Matrix(matrix);
		// A.eig().getV().print(10, 6);
		double[][] result = A.eig().getV().getArray();
		return result;
	}

	/**
	 * 寻找主成分
	 * 
	 * @param prinmaryArray
	 *            原始二维数据数组
	 * @param eigenvalue
	 *            特征值二维数组
	 * @param eigenVectors
	 *            特征向量二维数组
	 * @return principalMatrix 主成分矩阵
	 */
	public static Matrix getPrincipalComponent(double[][] primaryArray, double[][] eigenvalue,
			double[][] eigenVectors) {
		Matrix A = new Matrix(eigenVectors);// 定义一个特征向量矩阵
		double[][] tEigenVectors = A.transpose().getArray();// 特征向量转置
		Map<Integer, double[]> principalMap = new HashMap<Integer, double[]>();// key=主成分特征值，value=该特征值对应的特征向量
		TreeMap<Double, double[]> eigenMap = new TreeMap<Double, double[]>(Collections.reverseOrder());// key=特征值，value=对应的特征向量；初始化为翻转排序，使map按key值降序排列
		double total = 0;// 存储特征值总和
		int index = 0, n = eigenvalue.length;
		double[] eigenvalueArray = new double[n];// 把特征值矩阵对角线上的元素放到数组eigenvalueArray里
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j)
					eigenvalueArray[index] = eigenvalue[i][j];
			}
			index++;
		}

		for (int i = 0; i < tEigenVectors.length; i++) {
			double[] value = new double[tEigenVectors[0].length];
			value = tEigenVectors[i];
			eigenMap.put(eigenvalueArray[i], value);
		}

		// 求特征总和
		for (int i = 0; i < n; i++) {
			total += eigenvalueArray[i];
		}
		// 选出前几个主成分
		double temp = 0;
		@SuppressWarnings("unused")
		int principalComponentNum = 0;// 主成分数
		List<Double> plist = new ArrayList<Double>();// 主成分特征值
		for (double key : eigenMap.keySet()) {
			if (temp / total <= threshold) {
				temp += key;
				plist.add(key);
				principalComponentNum++;
			}
		}
		// System.out.println("\n" + "当前阈值: " + threshold);
		//System.out.println("取得的主成分数: " + principalComponentNum + "\n");

		// 往主成分map里输入数据
		for (int i = 0; i < plist.size(); i++) {
			if (eigenMap.containsKey(plist.get(i))) {
				principalMap.put(i, eigenMap.get(plist.get(i)));
			}
		}

		// 把map里的值存到二维数组里
		double[][] principalArray = new double[principalMap.size()][];
		// primaryArray[0][0]=0-primaryArray[0][0];
		Iterator<Entry<Integer, double[]>> it = principalMap.entrySet().iterator();
		for (int i = 0; it.hasNext(); i++) {
			principalArray[i] = it.next().getValue();
		}

		Matrix principalMatrix = new Matrix(principalArray);
		return principalMatrix;
	}

	/**
	 * 矩阵相乘
	 * 
	 * @param primary
	 *            原始二维数组
	 * 
	 * @param matrix
	 *            主成分矩阵
	 * 
	 * @return result 结果矩阵
	 */
	public static Matrix getResult(double[][] primary, Matrix matrix) {
		Matrix primaryMatrix = new Matrix(primary);
		Matrix result = primaryMatrix.times(matrix.transpose());
		return result;
	}

	public static double[][] pca1(double[][] data) {
		data = changeAverageToZero(data);
		// 协方差矩阵
		double[][] cov = getVarianceMatrix(data);
		// 特征值矩阵
		// System.out.println("特征值矩阵：");
		double[][] eigmatrix = getEigenvalueMatrix(cov);
		// 特征向量矩阵
		// System.out.println("特征向量矩阵：");
		double[][] eigvectorMatrix = getEigenVectorMatrix(cov);
		Matrix pmMatrix = getPrincipalComponent(data, eigmatrix, eigvectorMatrix);
		// double [][] ss=pmMatrix.getArray();
		// System.out.println("主成分矩阵：");
		double[][] newdata = getResult(data, pmMatrix).getArray();
		// System.out.println("新的特征：");
		return newdata;
	}

	/**
	 * PCA 降维
	 * 
	 * @param data
	 * @return 对于维数特别大，采用分割数据处理
	 */
	public static double[][] pca(double[][] data) {
		int length = data[0].length;

		double[][] newdata;
		if (length <= 1000) {
			newdata = pca1(data);
			return newdata;
		} else {
			int num = (int) Math.ceil((double) length / colsize);
			newdata = null;
			for (int i = 1; i <= num; i++) {
				double[][] split = splitData(data, i, colsize);
				double[][] subpca = pca1(split);
				newdata = mergeData(newdata, subpca);
			}
			return newdata;
		}
	}

	/**
	 * 分割数据
	 * 
	 * @param data
	 * @param i
	 * @param colsize
	 * @return
	 */
	public static double[][] splitData(double[][] data, int i, int colsize) {
		int row = data.length;
		int length = data[0].length;
		double[][] subdata;
		if ((i - 1) * colsize + colsize - 1 > length - 1) {
			subdata = new double[row][length - (i - 1) * colsize];
			for (int k = 0; k < row; k++) {
				for (int j = (i - 1) * colsize; j < length; j++) {
					int p = j - (i - 1) * colsize;
					subdata[k][p] = data[k][j];
				}
			}
		} else {
			subdata = new double[row][i * colsize - (i - 1) * colsize];
			for (int k = 0; k < row; k++) {
				for (int j = (i - 1) * colsize; j < i * colsize; j++) {
					int p = j - (i - 1) * colsize;
					subdata[k][p] = data[k][j];
				}
			}
		}
		return subdata;
	}
}
