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
 * �㷨����:
 * 1)��ԭʼ���ݰ������n��m�о���X
 * 2)�������Ļ�����ÿһά�����ݶ���ȥ��ά�ľ�ֵ��ʹÿһά�ľ�ֵ��Ϊ0
 * 3)���Э�������
 * 4)���Э������������ֵ����Ӧ����������
 * 5)��������������Ӧ������ֵ��С�������°������гɾ���ȡǰk����ɾ���p
 * 6)Y=PX ��Ϊ��ά��kά�������
 */
public class PCA {

	private static double threshold = 0.95;// false ����ֵ��ֵ
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
	 * �ϲ�����
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
	 * ת�ö�ά����
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
	 * ʹÿ�������ľ�ֵΪ0
	 * 
	 * @param primary
	 *            ԭʼ��ά�������
	 * @return averageArray ���Ļ���ľ���
	 */
	public static double[][] changeAverageToZero(double[][] primary) {
		// ��
		int n = primary.length;
		// ��
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
	 * ����Э�������
	 * 
	 * @param matrix
	 *            ���Ļ���ľ���
	 * @return result Э�������
	 */
	public static double[][] getVarianceMatrix(double[][] matrix) {
		int n = matrix.length;// ����
		int m = matrix[0].length;// ����
		double[][] result = new double[m][m];// Э�������
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
	 * ������ֵ����
	 * 
	 * @param matrix
	 *            Э�������
	 * @return result ����������ֵ��ά�������
	 */
	public static double[][] getEigenvalueMatrix(double[][] matrix) {
		Matrix A = new Matrix(matrix);
		// ������ֵ��ɵĶԽǾ���,eig()��ȡ����ֵ
		// A.eig().getD().print(10, 6);//��ӡ��ʽ���ո�С�����λ
		double[][] result = A.eig().getD().getArray();
		return result;
	}

	/**
	 * ��׼������������������
	 * 
	 * @param matrix
	 *            Э�������
	 * @return result ��׼����Ķ�ά�������
	 */
	public static double[][] getEigenVectorMatrix(double[][] matrix) {
		Matrix A = new Matrix(matrix);
		// A.eig().getV().print(10, 6);
		double[][] result = A.eig().getV().getArray();
		return result;
	}

	/**
	 * Ѱ�����ɷ�
	 * 
	 * @param prinmaryArray
	 *            ԭʼ��ά��������
	 * @param eigenvalue
	 *            ����ֵ��ά����
	 * @param eigenVectors
	 *            ����������ά����
	 * @return principalMatrix ���ɷ־���
	 */
	public static Matrix getPrincipalComponent(double[][] primaryArray, double[][] eigenvalue,
			double[][] eigenVectors) {
		Matrix A = new Matrix(eigenVectors);// ����һ��������������
		double[][] tEigenVectors = A.transpose().getArray();// ��������ת��
		Map<Integer, double[]> principalMap = new HashMap<Integer, double[]>();// key=���ɷ�����ֵ��value=������ֵ��Ӧ����������
		TreeMap<Double, double[]> eigenMap = new TreeMap<Double, double[]>(Collections.reverseOrder());// key=����ֵ��value=��Ӧ��������������ʼ��Ϊ��ת����ʹmap��keyֵ��������
		double total = 0;// �洢����ֵ�ܺ�
		int index = 0, n = eigenvalue.length;
		double[] eigenvalueArray = new double[n];// ������ֵ����Խ����ϵ�Ԫ�طŵ�����eigenvalueArray��
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

		// �������ܺ�
		for (int i = 0; i < n; i++) {
			total += eigenvalueArray[i];
		}
		// ѡ��ǰ�������ɷ�
		double temp = 0;
		@SuppressWarnings("unused")
		int principalComponentNum = 0;// ���ɷ���
		List<Double> plist = new ArrayList<Double>();// ���ɷ�����ֵ
		for (double key : eigenMap.keySet()) {
			if (temp / total <= threshold) {
				temp += key;
				plist.add(key);
				principalComponentNum++;
			}
		}
		// System.out.println("\n" + "��ǰ��ֵ: " + threshold);
		//System.out.println("ȡ�õ����ɷ���: " + principalComponentNum + "\n");

		// �����ɷ�map����������
		for (int i = 0; i < plist.size(); i++) {
			if (eigenMap.containsKey(plist.get(i))) {
				principalMap.put(i, eigenMap.get(plist.get(i)));
			}
		}

		// ��map���ֵ�浽��ά������
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
	 * �������
	 * 
	 * @param primary
	 *            ԭʼ��ά����
	 * 
	 * @param matrix
	 *            ���ɷ־���
	 * 
	 * @return result �������
	 */
	public static Matrix getResult(double[][] primary, Matrix matrix) {
		Matrix primaryMatrix = new Matrix(primary);
		Matrix result = primaryMatrix.times(matrix.transpose());
		return result;
	}

	public static double[][] pca1(double[][] data) {
		data = changeAverageToZero(data);
		// Э�������
		double[][] cov = getVarianceMatrix(data);
		// ����ֵ����
		// System.out.println("����ֵ����");
		double[][] eigmatrix = getEigenvalueMatrix(cov);
		// ������������
		// System.out.println("������������");
		double[][] eigvectorMatrix = getEigenVectorMatrix(cov);
		Matrix pmMatrix = getPrincipalComponent(data, eigmatrix, eigvectorMatrix);
		// double [][] ss=pmMatrix.getArray();
		// System.out.println("���ɷ־���");
		double[][] newdata = getResult(data, pmMatrix).getArray();
		// System.out.println("�µ�������");
		return newdata;
	}

	/**
	 * PCA ��ά
	 * 
	 * @param data
	 * @return ����ά���ر�󣬲��÷ָ����ݴ���
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
	 * �ָ�����
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
