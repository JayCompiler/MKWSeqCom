package seq.comp;

/**
 * Distance
 * @author zy
 *
 */
public class Distance {
	/*weighted Euclidean distance*/
	public static double euclidWeighted(double[] X,double[]Y,double[] weighted) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for ( int i=0;i<X.length;i++) {
			sum+=Math.pow((X[i]*weighted[i]-Y[i]*weighted[i]),2);
		}
		return Math.sqrt(sum);
	}
	/*weighted Manhattan distance*/
	public static double manhattanWeighted(double [] X,double[] Y,double[] weighted) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for (int i=0;i<X.length;i++) {
			sum+=Math.abs(X[i]-Y[i])*weighted[i];
		}
		return sum;
	}
	
	/*weighted Relative entropy*/
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
		return sum;
	}


	/**
	 * weighted Chebyshev distance
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double chebyshevWeighted(double[] X,double [] Y,double[]weight) {
		double max=-2147483647;
		for (int i = 0; i < Y.length; i++) {
			max=Math.max(max, Math.abs(X[i]*weight[i]-Y[i]*weight[i]));
		}
		return max;
	}
	
	
	/**
	 * weighted D2 similarity
	 * @param X k-mer count
	 * @param Y k-mer count
	 * @return
	 */
	public static double D2Weighted(double[] X,double [] Y,double [] weight) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i]*weight[i];
		}
		return 1/sum;
	}
	
	
	/**
	 * weighted D2S 
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
		return 1/sum;
	}
	/**
	 *weighted  D2* 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2starWeighted(double[] X,double [] Y,double[] weight) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i]*weight[i];
		}
		return 1/sum;
	}
	
	
	
	
	/*Euclidean distance*/
	public static double euclid(double[] X,double[]Y) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for ( int i=0;i<X.length;i++) {
			sum+=Math.pow((X[i]-Y[i]),2);
		}
		return Math.sqrt(sum);
	}
	/*Manhattan distance*/
	public static double manhattan(double [] X,double[] Y) {
		if(X.length!=Y.length)
			return -1.0;
		double sum=0.0;
		for (int i=0;i<X.length;i++) {
			sum+=Math.abs(X[i]-Y[i]);
		}
		return sum;
	}
	/*Relative entropy*/
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
		return sum;
	}
	/**
	 * Pearson correlation coefficient 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double pcc(double []X,double[] Y) {
		double pcc=Statistics.cov(X, Y)/(Statistics.sd(X)*Statistics.sd(Y));
		return 1/pcc;
	}
	/**
	 * Cosine distance
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
		return 1/cosine;
	}
	/**
	 * Chebyshev distance
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double chebyshev(double[] X,double [] Y) {
		double max=-2147483647;
		for (int i = 0; i < Y.length; i++) {
			max=Math.max(max, Math.abs(X[i]-Y[i]));
		}
		return max;
	}
	/**
	 * hao distance need to be smooth
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
		return sum;
	}
	/**
	 * D2 distance
	 * @param X k-mer count
	 * @param Y k-mer count
	 * @return
	 */
	public static double D2(double[] X,double [] Y) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i];
		}
		return 1/sum;
	}
	
	
	/**
	 * D2S 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2S(double[] X,double [] Y) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			double tmp =Math.sqrt((Math.pow(X[i], 2)+Math.pow(Y[i], 2)));
			sum+=(X[i]*Y[i])/tmp;
		}
		return 1/sum;
	}
	/**
	 * D2* 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double D2star(double[] X,double [] Y) {
		double sum=0.0;
		for (int i = 0; i < Y.length; i++) {
			sum+=X[i]*Y[i];
		}
		return 1/sum;
	}

}
