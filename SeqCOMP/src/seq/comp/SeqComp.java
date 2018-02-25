package seq.comp;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

public class SeqComp {
	private static final char[] dict = { 'A', 'T', 'C', 'G' };

	/**
	 * 提取序列特征
	 * 
	 * @param sequence
	 *            给定序列
	 * @param k
	 *            k-mer 中的k值
	 * @param allk_mer
	 *            所有可能的k_mer
	 * @param flag
	 * @return 提取到的k-mer t 特征 频率
	 */
	public static double[] seqFeaturefreq(String sequence, int k, String[] allk_mer,boolean flag) {
		String[] seqk_mer = getseqk_mer(sequence, k);
		/* 统计k_mer 频率 */
		HashMap<String, Double> xfreq = freq(seqk_mer, allk_mer,flag);
		double[] xf = new double[allk_mer.length];
		xf = xfreq2double(xfreq, allk_mer);
		//xf=normdata(xf);
		return xf;
	}
	/**
	 * 
	 * @param sequence 给定序列
	 * @param k k-mer 中的k值
	 * @param allk_mer   所有可能的k_mer
	 * @param flag  是否特征平滑
	 * @return
	 */
	public static double[] seqFeatureCount(String sequence, int k, String[] allk_mer,boolean flag) {
		String[] seqk_mer = getseqk_mer(sequence, k);
		/* 统计k_mer 频率 */
		HashMap<String, Double> xfreq = countkmer(seqk_mer, allk_mer,flag);
		double[] xf = new double[allk_mer.length];
		xf = xcount2double(xfreq, allk_mer);
		return xf;
	}

	/**
	 * json形式的特征值频率转换成double数组
	 * 
	 * @param xfreq
	 *            hashmap类型的特征值
	 * @param allk_mer
	 *            所有可能的k-mer
	 * @return 返回特征值double数组
	 */
	private static double[] xfreq2double(HashMap<String, Double> xfreq, String[] allk_mer) {
		double[] temp = new double[allk_mer.length];
		for (int i = 0; i < allk_mer.length; i++) {
			temp[i] = xfreq.get(allk_mer[i]);
		}
		return temp;
	}

	/**
	 * json形式的特征值频数转换成double数组
	 * 
	 * @param count
	 *            频数的 hashmap存储格式
	 * @param allk_mer
	 *            全部可能的k-mer
	 * @return
	 */
	public static double[] xcount2double(HashMap<String, Double> count, String[] allk_mer) {
		double[] temp = new double[allk_mer.length];
		for (int i = 0; i < allk_mer.length; i++) {
			temp[i] = count.get(allk_mer[i]);
		}
		return temp;
	}

	/**
	 * 
	 * @param k_mer
	 *            序列中提取到的 k-mer
	 * @param allk_mer
	 *            所有可能的k-mer
	 * @param flag 是否对数据平滑
	 * @return 统计序列k-mer频率 hashmap形式
	 */
	public static HashMap<String, Double> freq(String[] k_mer, String[] allk_mer, boolean flag) {
		double alpha = Math.pow(10, -8);
		double n1 = k_mer.length + allk_mer.length * alpha;
		double n2 = k_mer.length;
		HashMap<String, Double> kmerCount = countkmer(k_mer, allk_mer, flag);
		Iterator<Entry<String, Double>> iterator = kmerCount.entrySet().iterator();
		HashMap<String, Double> freq = new HashMap<>();
		while (iterator.hasNext()) {
			Map.Entry<String, Double> entry = (Map.Entry<String, Double>) iterator.next();
			if (flag) {
				freq.put(entry.getKey(), ((double) entry.getValue() / n1));
			} else {
				freq.put(entry.getKey(), ((double) entry.getValue() / n2));
			}
		}
		return freq;
	}

	/**
	 * 
	 * @param k_merSet
	 *            序列中提取到的k-mer
	 * @param allkmer
	 *            所有可能的 k-mer
	 * @param flag 
	 * 			是否对数据平滑
	 * @return 返回hashmap k_mer 频数 统计序列中k_mer的个数 ,平滑以后，即初始每个k-mer出现次数为alpha，不出现0的情况。
	 */
	public static HashMap<String, Double> countkmer(String[] k_merSet, String[] allkmer, boolean flag) {
		HashMap<String, Double> count = new HashMap<>();
		/* 初始化hashmap 平滑因子alpha */
		double alpha = Math.pow(10, -8);
		if (flag) {
			for (String k_mer : allkmer) 
				count.put(k_mer, alpha);// 设为alpha
		} else {
			for (String k_mer : allkmer) 
				count.put(k_mer, 0.0);// 若不平滑 设为0
		}
		for (String x : k_merSet) {
			if (count.get(x)==null) {
				System.out.println(x);
			}
			count.put(x, count.get(x) + 1.0);
		}
		return count;
	}

	/**
	 * @param dict
	 * @param k
	 * @return 将获得的k_mer转换成字符串数组
	 */
	public static String[] getK_mer(int k) {
		String b = new String();
		StringBuffer result = new StringBuffer();
		result = getAllK_mer(dict, k, b, result);
		return result.toString().split(" ");
	}

	/**
	 * 
	 * @param x
	 * @param k
	 *            k-mer中的k
	 * @return xk_mer 获得序列中的存在的k-mer
	 */
	public static String[] getseqk_mer(String X, int k) {
		int length = X.length();
		String[] xk_mer = new String[length - k + 1];
		for (int i = 0; i < length - k + 1; i++) {
			xk_mer[i] = X.substring(i, i + k);
		}
		return xk_mer;
	}

	/**
	 * 
	 * @param A
	 *            字符集
	 * @param k
	 *            k_mer中的k
	 * @param b
	 *            临时变量
	 * @param result
	 * @return result 获得 所有的k-mer O(4^k)
	 */
	public static StringBuffer getAllK_mer(char[] A, int k, String b, StringBuffer result) {
		if (k == 0) {
			result.append(b + " ");
		} else {
			for (int j = 0; j < A.length; j++) {
				b += A[j];
				getAllK_mer(A, k - 1, b, result);
				b = b.substring(0, b.length() - 1);
			}
		}
		return result;
	}
	
	public static double[] haoFeature(String sequence, int k) {
		/* 得到k-mer 和allk_mer */
		String[] k_Allmer = getK_mer(k);
		double[] k_1freq =seqFeaturefreq(sequence, k, k_Allmer, false);
		double[] prep=preProhao(sequence, k);
		double[] p=new double[(int)Math.pow(4, k)];
		for(int i=0;i<p.length;i++) {
			if(prep[i]==0) {
				p[i]=0;
			}else {
				p[i]=(k_1freq[i]-prep[i])/prep[i];
			}
		}
	
		return 	normdata(p);
	}

	/**
	 * 得到预测概率
	 * 
	 * @param sequence
	 *            基因序列
	 * @param k
	 * @return
	 */
	public static double[] preProhao(String sequence, int k) {
		/* 得到k-mer 和allk_mer */
		String[] k_Allmer = getK_mer( k);
		HashMap<String, Double> kfreq = new HashMap<>();
		kfreq = inithashmap(kfreq, k_Allmer);
		/* 得到(k-1)-mer,和all_k-1mer */
		String[] k_1mer = getseqk_mer(sequence, k - 1);
		String[] k_1Allmer = getK_mer( k - 1);
		HashMap<String, Double> k_1freq = freq(k_1mer, k_1Allmer,false);
		if (k < 1) {
			return null;
		}else if (k == 1) {
			return Markov.initProb(sequence);
		}else if (k == 2) {
			for (int i = 0; i < k_Allmer.length; i++) {
				double tmp = k_1freq.get(k_Allmer[i].substring(0, 1)) * k_1freq.get(k_Allmer[i].substring(1, 2));
				kfreq.put(k_Allmer[i], tmp);
			}
		} else {
			/* 得到(k-2)-mer和all_k-2mer */
			String[] k_2mer = getseqk_mer(sequence, k - 2);
			String[] k_2Allmer = getK_mer( k - 2);
			HashMap<String, Double> k_2freq = freq(k_2mer, k_2Allmer,true);
			for (int i = 0; i < k_Allmer.length; i++) {
				double tmp1 = k_1freq.get(k_Allmer[i].substring(0, k - 1)) * k_1freq.get(k_Allmer[i].substring(1, k));
				double tmp2 = tmp1 / k_2freq.get(k_Allmer[i].substring(1, k - 1));
				kfreq.put(k_Allmer[i], tmp2);
			}
		}
		int n = (int) Math.pow(4, k);
		double[] prep = new double[n];
		for (int i = 0; i < prep.length; i++) {
			String tmp = k_Allmer[i];
			prep[i] = kfreq.get(tmp);
		}
		return prep;
	}
		
	/**
	 * 仅对于D2S D2*距离 需要平滑数据
	 * @param sequence
	 * @param k
	 * @param allk_mer,
	 * @return
	 */
	public static double[] D2SCount(String sequence,int k,String[] allk_mer,int r) {
		int n=sequence.length()-k+1;
		double [] sc=seqFeatureCount(sequence, k, allk_mer, true);
		//double [] sc=seqFeaturefreq(sequence, k, allk_mer, true);
		//数据标准化
		double [] p=Markov.k_merMarkovProb(sequence, k,r);
		for(int i=0;i<sc.length;i++) {
			sc[i]=sc[i]-n*p[i];
		}
		return sc;
	}
	
	public static double[] D2starCount(String sequence,int k,String[] allk_mer,int r) {
		int n=sequence.length()-k+1;
		double [] sc=D2SCount(sequence, k, allk_mer, r);
		//sc=normdata(sc);
		//数据标准化
		double [] p=Markov.k_merMarkovProb(sequence, k,r);
		for(int i=0;i<sc.length;i++) {
			//double tmp =Math.sqrt(n*p[i]*(1-p[i]));
			double tmp =Math.sqrt(n*p[i]);
			if(tmp==0) {
				tmp=Math.pow(10, -8);
			}
			sc[i]=sc[i]/tmp;
		}
		return sc;
	}
	/**
	 * 多 k提取特征,采用zero-score
	 * @param sequence
	 * @param flag
	 * @return
	 */
	public static double[] seqAllFeaturefreq(String sequence,int startK,int kmax,boolean flag) {
		double[] feature=null;
//		double totalnum=0.0;
//		for(int k=startK;k<=kmax;k++) {
//			totalnum+=Math.pow(4, k);
//		}
		for(int k=startK;k<=kmax;k++) {
			String[] allk_mer=getK_mer(k);
			//double coe= Math.pow(4,k)/totalnum;
			//coe=Statistics.logarithm(2, 1/coe);
			double[] tmp=seqFeaturefreq(sequence, k, allk_mer, flag);
//			tmp=multiply(coe, tmp);
			tmp=normdata(tmp);
			feature=merge2double(tmp,feature);
		}
		//xf=normdata(xf);
		//改
		//feature=normdata(feature);
		return feature;
	}
	/**
	 * 全部频率，采用最大值最小值归一化。即没有负的特征
	 * @param sequence
	 * @param startK
	 * @param kmax
	 * @param flag
	 * @return
	 */
	public static double[] seqAllFeaturefreq2(String sequence,int startK,int kmax,boolean flag) {
		double[] feature=null;
//		double totalnum=0.0;
//		for(int k=startK;k<=kmax;k++) {
//			totalnum+=Math.pow(4, k);
//		}
		for(int k=startK;k<=kmax;k++) {
			String[] allk_mer=getK_mer(k);
			//double coe= Math.pow(4,k)/totalnum;
			//coe=Statistics.logarithm(2, 1/coe);
			double[] tmp=seqFeaturefreq(sequence, k, allk_mer, flag);
//			tmp=multiply(coe, tmp);
			tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//xf=normdata(xf);
		//改
		//feature=normdata(feature);
		return feature;
	}
	
	/**
	 * 多 k提取特征 D2
	 * @param sequence
	 * @param kmax
	 * @param flag
	 * @return
	 */
	public static double[] seqAllFeatureCount(String sequence,int startk, int kmax,boolean flag) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startk;k<=kmax;k++) {
			String[] allk_mer=getK_mer(k);
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			double[] tmp=seqFeatureCount(sequence, k, allk_mer, flag);
			//tmp=multiply(coe, tmp);
			//改1
			tmp=normdata(tmp);
			//tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//改
		//feature=normdata(feature);
		return feature;
	}
	/**
	 *  全部频率，采用最大值最小值归一化。即没有负的特征
	 * @param sequence
	 * @param startk
	 * @param kmax
	 * @param flag
	 * @return
	 */
	public static double[] seqAllFeatureCount2(String sequence,int startk, int kmax,boolean flag) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startk;k<=kmax;k++) {
			String[] allk_mer=getK_mer(k);
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			double[] tmp=seqFeatureCount(sequence, k, allk_mer, flag);
			//tmp=multiply(coe, tmp);
			//改1
			tmp=normdata2(tmp);
			//tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//改
		//feature=normdata(feature);
		return feature;
	}
	/**
	 * 多 k提取特征 D2S
	 * @param sequence
	 * @param kmax
	 * @param r
	 * @return
	 */
	public static double[] D2SAllCount(String sequence,int startK,int kmax,int r) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startK;k<=kmax;k++) {
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			String[] allk_mer=getK_mer(k);
			double[] tmp=D2SCount(sequence, k, allk_mer, r);
			//tmp=multiply(coe, tmp);
			//改 1
			tmp=normdata(tmp);
			//tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//改
		//feature=normdata(feature);
		return feature;
	}
	/**
	 * 采用最大值最小值归一化，不会出现负值
	 * @param sequence
	 * @param startK
	 * @param kmax
	 * @param r
	 * @return
	 */
	public static double[] D2SAllCount2(String sequence,int startK,int kmax,int r) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startK;k<=kmax;k++) {
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			String[] allk_mer=getK_mer(k);
			double[] tmp=D2SCount(sequence, k, allk_mer, r);
			//tmp=multiply(coe, tmp);
			//改 1
			tmp=normdata2(tmp);
			//tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//改
		//feature=normdata(feature);
		return feature;
	}
	/**
	 * 多 k提取特征
	 * @param sequence
	 * @param kmax
	 * @return
	 */
	public static double[] haoAllFeature(String sequence,int startK, int kmax) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startK;k<=kmax;k++) {
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			double[] tmp=haoFeature(sequence, k);
			//tmp=multiply(coe, tmp);
			tmp=normdata(tmp);
			feature=merge2double(tmp,feature);
		}
		return 	feature;
	}
	/**
	 * 最大值最小值归一化
	 * @param sequence
	 * @param startK
	 * @param kmax
	 * @return
	 */
	public static double[] haoAllFeature2(String sequence,int startK, int kmax) {
		double[] feature=null;
		//double totalnum=0.0;
		for(int k=startK;k<=kmax;k++) {
			//double coe=Math.pow(4,k)/totalnum;
			//coe=Statistics.logarithm(2, 1/coe);
			double[] tmp=haoFeature(sequence, k);
			//tmp=multiply(coe, tmp);
			tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		return 	feature;
	}
	/**
	 * @param sequence
	 * @param kmax
	 * @return
	 */
	public static double[] D2starAllCount(String sequence,int startK,int kmax,int r) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startK;k<=kmax;k++) {
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			String[] allk_mer=getK_mer(k);
			double[] tmp=D2starCount(sequence, k, allk_mer, r);
			//tmp=multiply(coe, tmp);
			//改 1
			tmp=normdata(tmp);
			//tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//改
		//feature=normdata(feature);
		return feature;
	}
	/**
	 * 采用最大值最小值归一化
	 * @param sequence
	 * @param startK
	 * @param kmax
	 * @param r
	 * @return
	 */
	public static double[] D2starAllCount2(String sequence,int startK,int kmax,int r) {
		double[] feature=null;
		double totalnum=0.0;
		for(int k=startK;k<=kmax;k++) {
			double coe=Math.pow(4,k)/totalnum;
			coe=Statistics.logarithm(2, 1/coe);
			String[] allk_mer=getK_mer(k);
			double[] tmp=D2starCount(sequence, k, allk_mer, r);
			//tmp=multiply(coe, tmp);
			//改 1
			tmp=normdata2(tmp);
			//tmp=normdata2(tmp);
			feature=merge2double(tmp,feature);
		}
		//改
		//feature=normdata(feature);
		return feature;
	}
	/*
	 * 数据乘上系数
	 */
	public static double[] multiply(double k,double[] data) {
		for(int i=0;i<data.length;i++) {
			data[i]=data[i]*k;
		}
		return data;
	}
	/**
	 * 合并数据
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[] merge2double(double[] X,double Y[]) {
		if(X==null&&Y!=null) {
			return Y;
		}else if(Y==null&&X!=null) {
			return X;
		}else if(X.length==0&&Y.length!=0) {
			return Y;
		}else if(Y.length==0&&X.length!=0) {
			return X;
		}else if(X.length!=0&&Y.length!=0) {
			int xlength=X.length;
			int ylength=Y.length;
			int zlength=xlength+ylength;
			int count=0;
			double [] Z=new double[zlength];
			for(int i=0;i<xlength;i++) {
				Z[count++]=X[i];
			}
			for(int j=0;j<ylength;j++) {
				Z[count++]=Y[j];
			}
			return Z;
		}else {
			return null;
		}
	}
	/**
	 * 标准化数据
	 * @param feature
	 * @return
	 */
	public static double[] normdata(double[] feature) {
		double mean=Statistics.mean(feature);
		double sd=Statistics.sd(feature);
		for(int i=0;i<feature.length;i++) {
			feature[i]=(feature[i]-mean)/sd;
		}
		return feature;
	}
	/**
	 * 最大值最小值归一化
	 * @param feature
	 * @return
	 */
	public static double[] normdata2(double[] feature) {
		double max=0.0;
		double min=Double.MAX_VALUE;
		for(int i=0;i<feature.length;i++) {
			if(feature[i]>max) {
				max=feature[i];
			}
			if(feature[i]<min) {
				min=feature[i];
			}
		}
		for(int i=0;i<feature.length;i++) {
			feature[i]=(feature[i]-min)/(max-min);
			if(feature[i]<0) {
				System.out.println("错错错错》》》》》》》》》》》》》》》");
			}
		}
		return feature;
	}
	/**
	 *  <string,double>hashmap 初始化
	 *  
	 */
	private static HashMap<String, Double> inithashmap(HashMap<String, Double> hashMap, String[] allkmer) {
		for (String string : allkmer) {
			hashMap.put(string, 0.0);
		}
		return hashMap;
	}
	/**
	 * 计算权重 基于最大离差方法
	 * @param freqs
	 * @return
	 */
	public static double[] getWeight(double[][] freqs) {
		int size=freqs[0].length;
		double[] weight=new double[size];
		double[] sumf=new double[size];
		for(int i=0;i<freqs[0].length;i++) {
			for(int j=0;j<freqs.length;j++) {
				for(int k=0;k<freqs.length;k++) {
					sumf[i]+=Math.abs(freqs[j][i]-freqs[k][i]);
				}
			}
		}
		double total=0.0;
		for(int i=0;i<size;i++) {
			total+=sumf[i];
		}
		for(int j=0;j<size;j++) {
			weight[j]=sumf[j]/total;
		}
		return weight;
	}
	
	
	/**
	 * 计算权重 基于熵值方法
	 * @param freqs
	 * @return
	 */
	public static double[] getWeight1(double[][] freqs) {
		for(int i=0;i<freqs.length;i++) {
			freqs[i]=normdata2(freqs[i]);
		}
		//行
		int row=freqs.length;
		//列
		int col=freqs[0].length;
		double K=1/Statistics.logarithm(Math.E, row);
		//System.out.println(K);
		//计算每列和
		double[] sum=new double[col];
		for (int i=0;i<col;i++) {
			for(int j=0;j<row;j++) {
				sum[i]=sum[i]+freqs[j][i];
			}
			if(sum[i]==0) {
				sum[i]=Math.pow(10, -10);
			}
		}
		double [][] contri=new double[row][col];
		//计算贡献度
		for(int i=0;i<col;i++) {
			for(int j=0;j<row;j++) {
				if(freqs[j][i]!=0) {
					contri[j][i]=freqs[j][i]/sum[i];
				}else {
					contri[j][i]=Math.pow(10, -10)/sum[i];
				}
			}
		}
		for(int i=0;i<col;i++) {
			for(int j=0;j<row;j++) {
				contri[j][i]=contri[j][i]*Statistics.logarithm(Math.E, contri[j][i]);
			}
		}

		//计算总贡献度
		double[] contribution=new double [col];
		for(int i=0;i<col;i++) {
			for(int j=0;j<row;j++) {
				contribution[i]+=-contri[j][i]*K;
			}
		}
		//计算一致性程度
		double [] consistency =new double[col];
		for(int i=0;i<col;i++) {
			consistency[i]=1-contribution[i];
		}
		//计算权重
		double [] weight=new double[col];
		double tmpSum=0.0;
		for(int i=0;i<col;i++) {
			tmpSum+=consistency[i];
		}
		for(int i=0;i<col;i++) {
			weight[i]=consistency[i]/tmpSum;
		}
		System.out.println();
		return weight;
	}
}
