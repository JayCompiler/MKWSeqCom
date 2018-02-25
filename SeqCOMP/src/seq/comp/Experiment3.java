package seq.comp;

import java.io.IOException;
import java.text.DecimalFormat;

public class Experiment3 {

	/**
	 * 测试数据 code代表提取特征的模式，有count ，freq， D2count haoFeature dis
	 * 表示距离模式：eur,cosine,ch,ma pcc，kld,hao,D2，D2S，D2star
	 * 
	 * eur,cosine,ch,ma ，pcc 不需要平滑 : flag=false code=0(freq); kld, 需要平滑 flag=true
	 * code=0 D2 flag =flase ,code=1 (count) D2S，需要平滑 flag=true code=2(D2count)
	 * r=0或1或2 hao特征 需要平滑 flag=true code =3 D2star需要平滑 flag=true code=4(D2starcount)
	 * r=0或1或2
	 */
	public static void testMultK() {
		// set parameter
		String dis = "D2star";
		int startK = 2;
		int kmax = 8;
		boolean flag = true;
		int code = 4;
		int r = 2;

		String[] data = ReadFile.readFileData2("src/e3/primates.txt", 27);
		String[] name = ReadFile.readFileDataName("src/e3/primates.txt", 27);

		for (int i = 0; i < name.length; i++) {
			if (name[i].length() > 10) {
				name[i] = name[i].substring(0, 10);
			} else {
				StringBuffer sBuffer = new StringBuffer();
				int c = 10 - name[i].length();
				if (c > 0) {
					for (int j = 0; j < c; j++) {
						sBuffer.append(" ");
					}
					name[i] = name[i] + sBuffer.toString();
				}
			}
		}

		double[][] freqs = new double[data.length][];
		for (int i = 0; i < data.length; i++) {
			freqs[i] = SeqComp.seqAllFeaturefreq(data[i], startK, kmax, flag);
		}
		double[] weight = SeqComp.getWeight(freqs);
		double[][] disMatrix = new double[data.length][data.length];

		// extract features
		double[][] features = new double[data.length][];
		if (code == 1) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.seqAllFeatureCount(data[i], startK, kmax, flag);
			}
		} else if (code == 2) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.D2SAllCount(data[i], startK, kmax, r);
			}
		} else if (code == 4) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.D2starAllCount(data[i], startK, kmax, r);
			}
		} else {
			System.out.println("输入code错误，请确认是1，2，4中一个！");
		}

		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data.length; j++) {
				disMatrix[i][j] = cptSim(features[i], features[j], weight, dis);
			}
		}
		// 标准化数据
		// 找最大值
		double max = -100;
		double min = 10000;
		for (int i = 0; i < disMatrix.length; i++) {
			for (int j = 0; j < disMatrix.length; j++) {
				if (disMatrix[i][j] > max) {
					max = disMatrix[i][j];
				}
				if (disMatrix[i][j] < min) {
					min = disMatrix[i][j];
				}
			}
		}
		double diff = max - min;
		for (int i = 0; i < disMatrix.length; i++) {
			for (int j = 0; j < disMatrix.length; j++) {
				disMatrix[i][j] = (disMatrix[i][j] - min) / diff;
			}
		}
		for (int i = 0; i < disMatrix.length; i++) {
			disMatrix[i][i] = 0.0;
		}
		DecimalFormat df = new DecimalFormat("0.000000");
		for (int i = 0; i < disMatrix.length; i++) {
			for (int j = 0; j < disMatrix.length; j++) {
				disMatrix[i][j] = Double.parseDouble(df.format(disMatrix[i][j]));
			}
		}

		//
		try {
			WriteFile.writeFile(disMatrix, name, dis + "-" + kmax + "-" + r);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * 测试单个k 数据 code代表提取特征的模式，有count ，freq， D2count haoFeature dis
	 * 表示距离模式：eur,cosine,ch,ma pcc，kld,hao,D2，D2S，D2star
	 * 
	 * eur,cosine,ch,ma ，pcc 不需要平滑 : flag=false code=0(freq); kld, 需要平滑 flag=true
	 * code=0 D2 flag =flase ,code=1 (count) D2S，需要平滑 flag=true code=2(D2count)
	 * r=0或1或2 hao特征 需要平滑 flag=true code =3 D2star需要平滑 flag=true code=4(D2starcount)
	 * r=0或1或2
	 */
	public static void testSigK() {
		// set parameter
		String dis = "D2star";
		int k =8;
		boolean flag = true;
		int code = 4;
		int r =2;
		String[] allk_mer = SeqComp.getK_mer(k);

		String[] data = ReadFile.readFileData2("src/e3/primates.txt", 27);
		String[] name = ReadFile.readFileDataName("src/e3/primates.txt", 27);
		for (int i = 0; i < name.length; i++) {
			if (name[i].length() > 10) {
				name[i] = name[i].substring(0, 10);
			} else {
				StringBuffer sBuffer = new StringBuffer();
				int c = 10 - name[i].length();
				if (c > 0) {
					for (int j = 0; j < c; j++) {
						sBuffer.append(" ");
					}
					name[i] = name[i] + sBuffer.toString();
				}
			}
		}
		double[][] disMatrix = new double[data.length][data.length];

		// extract features
		double[][] features = new double[data.length][];
		if (code == 1) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.seqFeatureCount(data[i], k, allk_mer, flag);
			}
		} else if (code == 2) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.D2SCount(data[i], k, allk_mer, r);
			}
		} else if (code == 4) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.D2starCount(data[i], k, allk_mer, r);
			}
		} else if (code == 0) {
			for (int i = 0; i < data.length; i++) {
				features[i] = SeqComp.seqFeaturefreq(data[i], k, allk_mer, flag);
			}
		} else {
			System.out.println("输入code错误，请确认是0,1，2，4中一个！");
		}

		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data.length; j++) {
				disMatrix[i][j] = cptSimSig(features[i], features[j], dis);
			}
		}
		// 标准化数据
		// 找最大值
		double max = -100;
		double min = 10000;
		for (int i = 0; i < disMatrix.length; i++) {
			for (int j = 0; j < disMatrix.length; j++) {
				if (disMatrix[i][j] > max) {
					max = disMatrix[i][j];
				}
				if (disMatrix[i][j] < min) {
					min = disMatrix[i][j];
				}
			}
		}
		double diff = max - min;
		for (int i = 0; i < disMatrix.length; i++) {
			for (int j = 0; j < disMatrix.length; j++) {
				disMatrix[i][j] = (disMatrix[i][j] - min) / diff;
			}
		}
		 for(int i=0;i<disMatrix.length;i++) {
		 disMatrix[i][i]=0.0;
		 }

		DecimalFormat df = new DecimalFormat("0.000000");
		for (int i = 0; i < disMatrix.length; i++) {
			for (int j = 0; j < disMatrix.length; j++) {
				disMatrix[i][j] = Double.parseDouble(df.format(disMatrix[i][j]));
			}
		}
		try {
			WriteFile.writeFile(disMatrix, name, dis + "-" + k + "-" + r);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param ds
	 * @param ds2
	 * @param dis
	 * @return code代表提取特征的模式，有count ，freq， D2count haoFeature dis
	 *         表示距离模式：eur,cosine,ch,ma pcc，kld,hao,D2，D2S，D2star
	 * 
	 *         eur,cosine,ch,ma ，pcc 不需要平滑 : flag=false code=0(freq); kld, 需要平滑
	 *         flag=true code=0 D2 flag =flase ,code=1 (count) D2S，需要平滑 flag=true
	 *         code=2(D2count) r=0或1或2 hao特征 需要平滑 flag=true code =3 D2star需要平滑
	 *         flag=true code=4(D2starcount) r=0或1或2
	 */

	private static double cptSimSig(double[] ds, double[] ds2, String dis) {
		double sim = 0.0;
		if (dis.equals("eur")) {
			sim = Distance.euclid(ds, ds2);
		} else if (dis.equals("ma")) {
			sim = Distance.manhattan(ds, ds2);
		} else if (dis.equals("cosine")) {
			sim = Distance.cosine(ds, ds2);
		} else if (dis.equals("ch")) {
			sim = Distance.chebyshev(ds, ds2);
		} else if (dis.equals("pcc")) {
			sim = Distance.pcc(ds, ds2);
		} else if (dis.equals("kld")) {
			sim = Distance.KLD(ds, ds2);
		} else if (dis.equals("D2")) {
			sim = Distance.D2(ds, ds2);
		} else if (dis.equals("D2S")) {
			sim = Distance.D2S(ds, ds2);
		} else if (dis.equals("D2star")) {
			sim = Distance.D2star(ds, ds2);
		} else {
			System.out.println("距离dis输入有误，请重新输入");
		}

		return sim;
	}

	private static double cptSim(double[] f1, double[] f2, double[] weight, String dis) {
		double sim = 0.0;
		if (dis.equals("D2")) {
			sim = Distance.D2Weighted(f1, f2, weight);
		} else if (dis.equals("D2S")) {
			sim = Distance.D2SWeighted(f1, f2, weight);
		} else if (dis.equals("D2star")) {
			sim = Distance.D2starWeighted(f1, f2, weight);
		}
		return sim;
	}

}
