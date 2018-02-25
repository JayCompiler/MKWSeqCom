package seq.comp;

public class Similarity {
	/*欧氏距离*/
	public static double euclidWeighted(double[] X,double[]Y,double[] weighted) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for ( int i=0;i<X.length;i++) {
			sum+=Math.pow((X[i]*weighted[i]-Y[i]*weighted[i]),2);
		}
		return 1/Math.sqrt(sum);
	}
	/*曼哈顿距离*/
	public static double manhattanWeighted(double [] X,double[] Y,double[] weighted) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for (int i=0;i<X.length;i++) {
			sum+=Math.abs(X[i]-Y[i])*weighted[i];
		}
		return 1/sum;
	}
	
	/*相对熵*/
	public static double KLDWeighted(double [] X,double[] Y,double[] weighted ){
		if(X.length!=Y.length)
			return -1.0;
		double k1=0.0;
		double k2=0.0;
		double sum=0.0;
		for (int i=0;i<X.length;i++) {
			k1+=X[i]*Statistics.logarithm(2, X[i]/Y[i]);
			k2+=Y[i]*Statistics.logarithm(2, Y[i]/X[i]);
		}
		sum=(k1+k2)/2;
		return 1/sum;
	}


	/**
	 * 切比雪夫距离
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double chebyshevWeighted(double[] X,double [] Y,double[]weight) {
		double max=-2147483647;
		for (int i = 0; i < Y.length; i++) {
			max=Math.max(max, Math.abs(X[i]*weight[i]-Y[i]*weight[i]));
		}
		return 1/max;
	}
	
	
	/**
	 * D2 相似度 
	 * @param X k-mer计数
	 * @param Y k-mer计数
	 * @return
	 */
	public static double D2Weighted(double[] X,double [] Y,double [] weight) {
		double sum=0.0;
//		X=SeqComp.normdata(X);
//		Y=SeqComp.normdata(Y);
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i]*weight[i];
		}
		return sum;
	}
	
	
	/**
	 * D2S 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2SWeighted(double[] X,double [] Y,double [] weight) {
		double sum=0.0;
//		X=SeqComp.normdata(X);
//		Y=SeqComp.normdata(Y);
		for (int i = 0; i < Y.length; i++) {
			double tmp =Math.sqrt((Math.pow(X[i], 2)+Math.pow(Y[i], 2)));
			sum+=(X[i]*Y[i]*weight[i])/tmp;
		}
		//sum=Math.abs(sum);
		return sum;
	}
	/**
	 * D2* 距离 越大越相关所以与距离保持一致取个倒数
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2starWeighted(double[] X,double [] Y,double[] weight) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i]*weight[i];
		}
		return sum;
	}
	
	
	
	
	/*欧氏距离*/
	public static double euclid(double[] X,double[]Y) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for ( int i=0;i<X.length;i++) {
			sum+=Math.pow((X[i]-Y[i]),2);
		}
		return 1/Math.sqrt(sum);
	}
	/*曼哈顿距离*/
	public static double manhattan(double [] X,double[] Y) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for (int i=0;i<X.length;i++) {
			sum+=Math.abs(X[i]-Y[i]);
		}
		return 1/sum;
	}
	/*相对熵*/
	public static double KLD(double [] X,double[] Y) {
		if(X.length!=Y.length)
			return -1.0;
		double k1=0.0;
		double k2=0.0;
		double sum=0.0;
		for (int i=0;i<X.length;i++) {
			k1+=X[i]*Statistics.logarithm(2, X[i]/Y[i]);
			k2+=Y[i]*Statistics.logarithm(2, Y[i]/X[i]);
		}
		sum=(k1+k2)/2;
		return 1/sum;
	}
	/**
	 * 皮尔僧相关系数 ,越大越相关所以与距离(保持一致取个倒数)
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double pcc(double []X,double[] Y) {
		double pcc=Statistics.cov(X, Y)/(Statistics.sd(X)*Statistics.sd(Y));
		return pcc;
	}
	/**
	 * 余弦相似度 越大越相关所以与距离(保持一致取个倒数)
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double cosine(double [] X,double [] Y) {
		double sum=0.0;
		double lx=0.0;
		double ly=0.0;
		double cosine=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i];
			lx+=Math.pow(X[i],2);
			ly+=Math.pow(Y[i],2);
		}
		cosine=sum/(Math.sqrt(lx)*Math.sqrt(ly));
		return cosine;
	}
	/**
	 * 切比雪夫距离
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double chebyshev(double[] X,double [] Y) {
		double max=-2147483647;
		for (int i = 0; i < Y.length; i++) {
			max=Math.max(max, Math.abs(X[i]-Y[i]));
		}
		return 1/max;
	}
	/**
	 * hao式距离 需要平滑
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double hao(double[] X,double [] Y) {
		double sum=0.0;
		double tp1=0.0;
		double tp2=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=(X[i]*Y[i]);
			tp1+=Math.pow(X[i], 2);
			tp2+=Math.pow(Y[i], 2);
		}
		if(Math.sqrt((tp1*tp2))==0) {
			sum=0;
		}else {
			sum=0.5*(1-sum/Math.sqrt((tp1*tp2)));
		}
		return 1/sum;
	}
	/**
	 * D2 相似度 
	 * @param X k-mer计数
	 * @param Y k-mer计数
	 * @return
	 */
	public static double D2(double[] X,double [] Y) {
		double sum=0.0;
//		X=SeqComp.normdata(X);
//		Y=SeqComp.normdata(Y);
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i];
		}
		return sum;
	}
	
	
	/**
	 * D2S 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2S(double[] X,double [] Y) {
		double sum=0.0;
//		X=SeqComp.normdata(X);
//		Y=SeqComp.normdata(Y);
		for (int i = 0; i < Y.length; i++) {
			double tmp =Math.sqrt((Math.pow(X[i], 2)+Math.pow(Y[i], 2)));
			sum+=(X[i]*Y[i])/tmp;
		}
	   // Math.abs(sum);
		return sum;
	}
	/**
	 * D2* 相似度
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2star(double[] X,double [] Y) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i];
		}
		return sum;
	}

}
