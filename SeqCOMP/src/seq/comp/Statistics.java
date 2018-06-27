package seq.comp;
/* 统计学的一些相关计算*/
public class Statistics {
	
	public static void main(String[] args) {
		double [] X= {5.0, 3.0, 2.5};
		double [] Y= {4.0, 3.0, 2.0};
		System.out.println(mean(X));
		System.out.println(var(X));
		System.out.println(sd(X));
		System.out.println(cov(X, Y));
	}
	/**
	 * 最大值求取
	 * @param X
	 * @return
	 */
	public static double max(double[] X) {
		double max=Double.MIN_VALUE;
		for(int i=0;i<X.length;i++) {
			if(X[i]>max) {
				max=X[i];
			}
		}
		return max;
	}
	/**
	 * @param X
	 * @return
	 * 计算平均数
	 */
	public static double mean(double[] X) {
		double sum=0.0;
		for (double d : X) {
			sum+=d;
		}
		return sum/X.length;
	}
	/**
	 * 计算无偏方差
	 * @param X
	 * @return
	 */
	public static double var(double[] X) {
		double sum=0.0;
		double mean=mean(X);
		for (double d : X) {
			sum+=Math.pow((d-mean),2);
		}
		return sum/(X.length-1);
	}
	/**
	 * 计算无偏标准差
	 * @param X
	 * @return
	 */
	public static double sd(double[] X) {
		double sd=Math.sqrt(var(X));
		return sd;	
	}
	/**
	 * 计算协方差
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double cov(double[] X,double [] Y) {
		double sum=0.0;
		double meanX=mean(X);
		double meanY=mean(Y);
		for (int i=0;i<X.length;i++) {
			sum+=(X[i]-meanX)*(Y[i]-meanY);
		}
		return sum/(X.length-1);
	}
	/**
	 * 计算对数
	 * @param base
	 * @param value
	 * @return
	 */
	public static double logarithm(double base,double value) {
		return Math.log(value)/Math.log(base);
	}

}
