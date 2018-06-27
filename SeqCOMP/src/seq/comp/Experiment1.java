package seq.comp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class Experiment1 {

	/**
	 * no pca ,test difference k,r
	 * 
	 * @param k
	 * @param r
	 * @return
	 */
	public static double[] testWeight(int k, int r, boolean flag, int code, String dis) {
		/* get 4^k k-mer array */
		String[] allk_mer = SeqComp.getK_mer(k);
		// 读取目录
		ArrayList<String> filename = ReadFile.readAllFile("src/dataset/");
		String querySequence = ReadFile.readFileData("src/dataset/HSLIPAS.fasta");
		//
		querySequence = ReadFile.standData(querySequence);
		String[] other = new String[39];
		double[][] otherFeatures = new double[39][(int) Math.pow(4, k)];
		double[][] otherfreqs = new double[39][];
		/**
		 * 相似度 结果result
		 */
		HashMap<String, Double> result = new HashMap<>();
		HashMap<String, String> map = new HashMap<>();
		for (int i = 0; i < filename.size(); i++) {
			String tmp1 = filename.get(i).split("\\.")[0].trim();
			if (tmp1.equals("HSLIPAS")) {

			} else {
				String tmp = "src/dataset/" + filename.get(i);
				String sequence = ReadFile.readFileData(tmp);
				map.put(filename.get(i).split("\\.")[0], sequence);
			}
		}
		// 获取查询序列的特征向量
		double[] queryFeature;
		double[] queryfreq = SeqComp.seqFeaturefreq(querySequence, k, allk_mer, flag);
		if (code == 0) {
			queryFeature = SeqComp.seqFeaturefreq(querySequence, k, allk_mer, flag);
		} else if (code == 1) {
			queryFeature = SeqComp.seqFeatureCount(querySequence, k, allk_mer, flag);
		}
		// D2,D2S,D2star 需要计算词频 D2S,D2*需要平滑
		else if (code == 2) {
			queryFeature = SeqComp.D2SCount(querySequence, k, allk_mer, r);
		} else if (code == 3) {
			queryFeature = SeqComp.haoFeature(querySequence, k);
		} else if (code == 4) {
			queryFeature = SeqComp.D2starCount(querySequence, k, allk_mer, r);
		} else {
			queryFeature = null;
			System.out.println("查询序列错误");
		}

		/**
		 * 获取其余序列的特征向量
		 */
		Iterator<Map.Entry<String, String>> it = map.entrySet().iterator();
		int count = 0;
		String[] others = new String[39];
		ArrayList<GenePair> genePairs = new ArrayList<>();
		GenePair[] genePair = new GenePair[39];

		while (it.hasNext()) {
			Map.Entry<String, String> entry = it.next();
			boolean pos = false;
			if (entry.getKey().indexOf("+") != -1) {
				pos = true;
			} else {
				pos = false;
			}

			String seq = ReadFile.standData(entry.getValue());
			if (code == 0) {
				// System.out.println(code);
				otherFeatures[count] = SeqComp.seqFeaturefreq(seq, k, allk_mer, flag);
			} else if (code == 1) {
				otherFeatures[count] = SeqComp.seqFeatureCount(seq, k, allk_mer, flag);
			} else if (code == 2) {
				otherFeatures[count] = SeqComp.D2SCount(seq, k, allk_mer, r);
			} else if (code == 3) {
				otherFeatures[count] = SeqComp.haoFeature(seq, k);
			} else if (code == 4) {
				otherFeatures[count] = SeqComp.D2starCount(seq, k, allk_mer, r);
			} else {
				// System.out.println("这是"+code);
				otherFeatures[count] = null;
				System.out.println("查询序列错误");
			}
			otherfreqs[count] = SeqComp.seqFeaturefreq(seq, k, allk_mer, flag);
			others[count] = seq;
			other[count] = entry.getKey();
			int length = querySequence.length() + seq.length();
			// 构造基因对
			genePair[count] = new GenePair();
			genePair[count].setName(entry.getKey());
			genePair[count].setLength(length);
			genePair[count].setPos(pos);
			genePair[count].setSequence1(querySequence);
			genePair[count].setSequence2(seq);
			count++;
		}

		double[][] tmpFeature = new double[40][];
		tmpFeature[0] = queryFeature;
		double[][] freqs = new double[40][];
		freqs[0] = queryfreq;
		for (int i = 1; i < tmpFeature.length; i++) {
			tmpFeature[i] = otherFeatures[i - 1];
			freqs[i] = otherfreqs[i - 1];
		}
		// 使用freqs 计算权重
		double[] weight = SeqComp.getWeight(freqs);
		// double[][] pcaFeature = tmpFeature;
		double maxweight = Statistics.max(weight);
		System.out.println("max WEight k="+k+":"+maxweight);
		// 计算相似度
		for (int i = 0; i < otherFeatures.length; i++) {
			double sim = 0.0;
			if (dis.equals("eur")) {
				sim = Similarity.euclid(queryFeature, otherFeatures[i]);
			} else if (dis.equals("cosine")) {
				sim = Similarity.cosine(queryFeature, otherFeatures[i]);
			} else if (dis.equals("ch")) {
				sim = Similarity.chebyshev(queryFeature, otherFeatures[i]);
			} else if (dis.equals("ma")) {
				sim = Similarity.manhattan(queryFeature, otherFeatures[i]);
			} else if (dis.equals("pcc")) {
				sim = Similarity.pcc(queryFeature, otherFeatures[i]);
			} else if (dis.equals("kld")) {
				sim = Similarity.KLD(queryFeature, otherFeatures[i]);
			} else if (dis.equals("hao")) {
				sim = Similarity.hao(queryFeature, otherFeatures[i]);
			} else if (dis.equals("D2")) {
				sim = Similarity.D2Weighted(queryFeature, otherFeatures[i], weight);
			} else if (dis.equals("D2S")) {
				sim = Similarity.D2SWeighted(queryFeature, otherFeatures[i], weight);
			} else if (dis.equals("D2star")) {
				sim = Similarity.D2starWeighted(queryFeature, otherFeatures[i], weight);
			} else {
				sim = 0.0;
				System.out.println("dis输入有误");
			}
			genePair[i].setTotal_sim(sim);
			genePairs.add(genePair[i]);
			result.put(other[i], sim);
		}

		// 将map.entrySet()转换成list
		List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(result.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String, Double>>() {
			// 降序排序
			@Override
			public int compare(Entry<String, Double> o1, Entry<String, Double> o2) {
				// return o1.getValue().compareTo(o2.getValue());
				return o2.getValue().compareTo(o1.getValue()); // 从大到小排列，相似度越大，越相似
			}
		});

		int count1 = 0;
		String[] resultS = new String[39];

		for (Map.Entry<String, Double> mapping : list) {
			String tmp = count1 + ":" + mapping.getKey() + ":" + mapping.getValue();
			resultS[count1] = tmp;
			count1++;
			// System.out.println(tmp);
		}
		double[] acc = accuracy(resultS);
		double[] rs = new double[5];
		rs[0] = acc[0];
		rs[1] = acc[1];
		rs[2] = acc[2];
		rs[3] = acc[3];
		rs[4] = calROCRes(genePairs);
		double[][] coodinate = calROCCood(genePairs);
		try {
			WriteFile.writeCood(coodinate, "Sigweight" + dis + "-" + k + "-" + r);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return rs;
	}

	/**
	 * 带权重的 多k
	 * 
	 * @param flag
	 * @param kmax
	 * @param code
	 * @param dis
	 * @param r
	 * @return
	 */
	public static double[] testAllkWeight(boolean flag, int startk, int kmax, int code, String dis, int r) {
		// 读取目录
		ArrayList<String> filename = ReadFile.readAllFile("src/dataset/");
		String querySequence = ReadFile.readFileData("src/dataset/HSLIPAS.fasta");
		//
		querySequence = ReadFile.standData(querySequence);
		String[] other = new String[39];
		double[][] otherFeatures = new double[39][];
		double[][] otherfreqs = new double[39][];
		/**
		 * 相似度 结果result
		 */
		HashMap<String, Double> result = new HashMap<>();
		HashMap<String, String> map = new HashMap<>();
		for (int i = 0; i < filename.size(); i++) {
			String tmp1 = filename.get(i).split("\\.")[0].trim();
			if (tmp1.equals("HSLIPAS")) {

			} else {
				String tmp = "src/dataset/" + filename.get(i);
				String sequence = ReadFile.readFileData(tmp);
				map.put(filename.get(i).split("\\.")[0], sequence);
			}
		}
		// 获取查询序列的特征向量
		double[] queryFeature;
		// 查询序列所有频率
		double[] queryfreq = SeqComp.seqAllFeatureCount(querySequence, startk, kmax, flag);
		if (code == 0) {
			queryFeature = SeqComp.seqAllFeaturefreq(querySequence, startk, kmax, flag);
		} else if (code == 1) {
			queryFeature = SeqComp.seqAllFeatureCount(querySequence, startk, kmax, flag);
		} else if (code == 2) {
			queryFeature = SeqComp.D2SAllCount(querySequence, startk, kmax, r);
		} else if (code == 3) {
			queryFeature = SeqComp.haoAllFeature(querySequence, startk, kmax);
		} else if (code == 4) {
			queryFeature = SeqComp.D2starAllCount(querySequence, startk, kmax, r);
		} else {
			queryFeature = null;
			System.out.println("查询序列错误");
		}

		//
		/**
		 * 获取其余序列的特征向量
		 */
		Iterator<Map.Entry<String, String>> it = map.entrySet().iterator();
		int count = 0;
		String[] others = new String[39];

		ArrayList<GenePair> genePairs = new ArrayList<>();

		GenePair[] genePair = new GenePair[39];

		while (it.hasNext()) {
			Map.Entry<String, String> entry = it.next();
			boolean pos = false;
			if (entry.getKey().indexOf("+") != -1) {
				pos = true;
			} else {
				pos = false;
			}

			String seq = ReadFile.standData(entry.getValue());
			if (code == 0) {
				otherFeatures[count] = SeqComp.seqAllFeaturefreq(seq, startk, kmax, flag);
			} else if (code == 1) {
				otherFeatures[count] = SeqComp.seqAllFeatureCount(seq, startk, kmax, flag);
			} else if (code == 2) {
				otherFeatures[count] = SeqComp.D2SAllCount(seq, startk, kmax, r);
			} else if (code == 3) {
				otherFeatures[count] = SeqComp.haoAllFeature(seq, startk, kmax);
			} else if (code == 4) {
				otherFeatures[count] = SeqComp.D2starAllCount(seq, startk, kmax, r);
			} else {
				otherFeatures[count] = null;
				System.out.println("查询序列错误");
			}
			otherfreqs[count] = SeqComp.seqAllFeaturefreq(seq, startk, kmax, flag);
			others[count] = seq;
			other[count] = entry.getKey();
			int length = querySequence.length() + seq.length();
			// 构造基因对
			genePair[count] = new GenePair();
			genePair[count].setName(entry.getKey());
			genePair[count].setLength(length);
			genePair[count].setPos(pos);
			genePair[count].setSequence1(querySequence);
			genePair[count].setSequence2(seq);
			count++;
		}
		// 加上权重
		double[][] tmpFeature = new double[40][];
		tmpFeature[0] = queryFeature;
		double[][] freqs = new double[40][];
		freqs[0] = queryfreq;
		for (int i = 1; i < tmpFeature.length; i++) {
			tmpFeature[i] = otherFeatures[i - 1];
			freqs[i] = otherfreqs[i - 1];
		}
		// 使用freqs 计算权重
		double[] weight = SeqComp.getWeight(freqs);
		int index=0;
		double[] w2=new double[(int)Math.pow(4, 2)];
		for(;index<(int)Math.pow(4, 2);index++) {
			w2[index]=weight[index];
		}
		double max2=Statistics.max(w2);
		double mean2=Statistics.mean(w2);
		double[] w3=new double[(int)Math.pow(4, 3)];
		for(;(index-(int)Math.pow(4, 2))<(int)Math.pow(4, 3);index++) {
			w3[index-(int)Math.pow(4, 2)]=weight[index];
		}
		double max3=Statistics.max(w3);
		double mean3=Statistics.mean(w3);
		double[] w4=new double[(int)Math.pow(4, 4)];
		for(;(index-(int)Math.pow(4, 3))<(int)Math.pow(4, 4);index++) {
			w4[index-(int)Math.pow(4, 3)]=weight[index];
		}
		double mean4=Statistics.mean(w4);
		double max4=Statistics.max(w4);
		double[] w5=new double[(int)Math.pow(4, 5)];
		for(;(index-(int)Math.pow(4, 4))<(int)Math.pow(4, 5);index++) {
			w5[index-(int)Math.pow(4, 4)]=weight[index];
		}
		double max5=Statistics.max(w5);
		double mean5=Statistics.mean(w5);
		double[] w6=new double[(int)Math.pow(4, 6)];
		for(;(index-(int)Math.pow(4, 5))<(int)Math.pow(4, 6);index++) {
			w6[index-(int)Math.pow(4, 5)]=weight[index];
		}
		double mean6=Statistics.mean(w6);
		double max6=Statistics.max(w6);
		double[] w7=new double[(int)Math.pow(4, 7)];
		for(;(index-(int)Math.pow(4, 6))<(int)Math.pow(4, 7);index++) {
			w7[index-(int)Math.pow(4, 6)]=weight[index];
		}
		double max7=Statistics.max(w7);
		double mean7=Statistics.mean(w7);
		double[] w8=new double[(int)Math.pow(4, 8)];
		for(;(index-(int)Math.pow(4, 7))<(int)Math.pow(4,8);index++) {
			w8[index-(int)Math.pow(4, 7)]=weight[index];
		}
		double mean8=Statistics.mean(w8);
		double max8=Statistics.max(w8);
		System.out.println("max2:"+max2+" max3:"+max3
				+" max4:"+max4+" max5:"+max5+
				" max6:"+max6+" max7:"+max7+" max8:"+max8);
		System.out.println("m2:"+mean2+" m3:"+mean3
				+" m4:"+mean4+" m5:"+mean5+
				" m6:"+mean6+" m7:"+mean7+" m8:"+mean8);
		double[][] pcaFeature = tmpFeature;

		// 计算相似度
		for (int i = 0; i < pcaFeature.length - 1; i++) {
			double sim = 0.0;
			if (dis.equals("eur")) {
				// sim = Similarity.euclidWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.euclid(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("cosine")) {
				sim = Similarity.cosine(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("ch")) {
				// sim = Similarity.chebyshevWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.chebyshev(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("ma")) {
				// sim = Similarity.manhattanWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.manhattan(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("pcc")) {
				sim = Similarity.pcc(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("kld")) {
				// sim = Similarity.KLDWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.KLD(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("hao")) {
				sim = Similarity.hao(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("D2")) {
				sim = Similarity.D2Weighted(pcaFeature[0], pcaFeature[i + 1], weight);
			} else if (dis.equals("D2S")) {
				sim = Similarity.D2SWeighted(pcaFeature[0], pcaFeature[i + 1], weight);
			} else if (dis.equals("D2star")) {
				sim = Similarity.D2starWeighted(queryFeature, otherFeatures[i], weight);
			} else {
				sim = 0.0;
				System.out.println("dis输入有误");
			}
			genePair[i].setTotal_sim(sim);
			genePairs.add(genePair[i]);
			result.put(other[i], sim);
		}

		// 将map.entrySet()转换成list
		List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(result.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String, Double>>() {
			// 降序排序
			@Override
			public int compare(Entry<String, Double> o1, Entry<String, Double> o2) {
				// return o1.getValue().compareTo(o2.getValue());
				return o2.getValue().compareTo(o1.getValue()); // 从大到小排列，相似度越大，越相似
			}
		});

		int count1 = 0;
		String[] resultS = new String[39];

		for (Map.Entry<String, Double> mapping : list) {
			String tmp = count1 + ":" + mapping.getKey() + ":" + mapping.getValue();
			resultS[count1] = tmp;
			count1++;
			// System.out.println(tmp);
		}
		double[] acc = accuracy(resultS);
		double[] rs = new double[5];
		rs[0] = acc[0];
		rs[1] = acc[1];
		rs[2] = acc[2];
		rs[3] = acc[3];
		rs[4] = calROCRes(genePairs);
		if (code == 1 || code == 2 || code == 4) {
			double[][] coodinate = calROCCood(genePairs);
			try {
				WriteFile.writeCood(coodinate, "MUlweight" + dis + "-" + kmax + "-" + r);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("the  " + dis + " multik is:");
		return rs;
	}

	/**
	 * 打印 hashmap
	 * 
	 * @param a
	 */
	public static void printhashmap(HashMap<String, Double> a) {
		Iterator<Entry<String, Double>> iterator = a.entrySet().iterator();
		while (iterator.hasNext()) {
			Map.Entry<String, Double> entry = (Entry<String, Double>) iterator.next();
			System.out.println(entry.getKey() + ":" + entry.getValue());
		}
	}

	/**
	 * 统计TP,FP,TN,FN
	 * 
	 * @param resultS
	 * @return
	 */
	public static double[] accuracy(String[] resultS) {
		double[] acc = new double[4];
		double TP = 0.0; // 预测相关正确
		double FP = 0.0; // 预测相关错误
		double TN = 0.0; // 预测不相关正确
		double FN = 0.0; // 预测不相关错误
		for (int i = 0; i < 20; i++) {
			if (resultS[i].indexOf("+:") != -1) {
				TP = TP + 1;
			}
		}
		FP = 20 - TP;
		FN = FP;
		TN = 19 - FN;
		acc[0] = TP;
		acc[1] = FP;
		acc[2] = TN;
		acc[3] = FN;
		return acc;
	}

	/**
	 * 不进行pca处理 测试不同的k,r
	 * 
	 * @param k
	 * @param r
	 * @return
	 */
	public static double[] test(int k, int r, boolean flag, int code, String dis) {
		/* 获得4^k个k-mer数组 */
		String[] allk_mer = SeqComp.getK_mer(k);
		// 读取目录
		ArrayList<String> filename = ReadFile.readAllFile("src/dataset/");
		String querySequence = ReadFile.readFileData("src/dataset/HSLIPAS.fasta");
		//
		querySequence = ReadFile.standData(querySequence);
		String[] other = new String[39];
		double[][] otherFeatures = new double[39][(int) Math.pow(4, k)];
		/**
		 * 相似度 结果result
		 */
		HashMap<String, Double> result = new HashMap<>();
		HashMap<String, String> map = new HashMap<>();
		for (int i = 0; i < filename.size(); i++) {
			String tmp1 = filename.get(i).split("\\.")[0].trim();
			if (tmp1.equals("HSLIPAS")) {

			} else {
				String tmp = "src/dataset/" + filename.get(i);
				String sequence = ReadFile.readFileData(tmp);
				map.put(filename.get(i).split("\\.")[0], sequence);
			}
		}
		// 获取查询序列的特征向量
		double[] queryFeature;
		if (code == 0) {
			queryFeature = SeqComp.seqFeaturefreq(querySequence, k, allk_mer, flag);
		} else if (code == 1) {
			queryFeature = SeqComp.seqFeatureCount(querySequence, k, allk_mer, flag);
		}
		// D2,D2S,D2star 需要计算词频 D2S,D2*需要平滑
		else if (code == 2) {
			queryFeature = SeqComp.D2SCount(querySequence, k, allk_mer, r);
		} else if (code == 3) {
			queryFeature = SeqComp.haoFeature(querySequence, k);
		} else if (code == 4) {
			queryFeature = SeqComp.D2starCount(querySequence, k, allk_mer, r);
		} else {
			queryFeature = null;
			System.out.println("查询序列错误");
		}

		/**
		 * 获取其余序列的特征向量
		 */
		Iterator<Map.Entry<String, String>> it = map.entrySet().iterator();
		int count = 0;
		String[] others = new String[39];

		ArrayList<GenePair> genePairs = new ArrayList<>();

		GenePair[] genePair = new GenePair[39];

		while (it.hasNext()) {
			Map.Entry<String, String> entry = it.next();
			boolean pos = false;
			if (entry.getKey().indexOf("+") != -1) {
				pos = true;
			} else {
				pos = false;
			}

			String seq = ReadFile.standData(entry.getValue());
			if (code == 0) {
				// System.out.println(code);
				otherFeatures[count] = SeqComp.seqFeaturefreq(seq, k, allk_mer, flag);
			} else if (code == 1) {
				otherFeatures[count] = SeqComp.seqFeatureCount(seq, k, allk_mer, flag);
			} else if (code == 2) {
				otherFeatures[count] = SeqComp.D2SCount(seq, k, allk_mer, r);
			} else if (code == 3) {
				otherFeatures[count] = SeqComp.haoFeature(seq, k);
			} else if (code == 4) {
				otherFeatures[count] = SeqComp.D2starCount(seq, k, allk_mer, r);
			} else {
				// System.out.println("这是"+code);
				otherFeatures[count] = null;
				System.out.println("查询序列错误");
			}
			others[count] = seq;
			other[count] = entry.getKey();
			int length = querySequence.length() + seq.length();
			// 构造基因对
			genePair[count] = new GenePair();
			genePair[count].setName(entry.getKey());
			genePair[count].setLength(length);
			genePair[count].setPos(pos);
			genePair[count].setSequence1(querySequence);
			genePair[count].setSequence2(seq);
			count++;
		}
		// 计算相似度
		for (int i = 0; i < otherFeatures.length; i++) {
			double sim = 0.0;
			if (dis.equals("eur")) {
				sim = Similarity.euclid(queryFeature, otherFeatures[i]);
			} else if (dis.equals("cosine")) {
				sim = Similarity.cosine(queryFeature, otherFeatures[i]);
			} else if (dis.equals("ch")) {
				sim = Similarity.chebyshev(queryFeature, otherFeatures[i]);
			} else if (dis.equals("ma")) {
				sim = Similarity.manhattan(queryFeature, otherFeatures[i]);
			} else if (dis.equals("pcc")) {
				sim = Similarity.pcc(queryFeature, otherFeatures[i]);
			} else if (dis.equals("kld")) {
				sim = Similarity.KLD(queryFeature, otherFeatures[i]);
			} else if (dis.equals("hao")) {
				sim = Similarity.hao(queryFeature, otherFeatures[i]);
			} else if (dis.equals("D2")) {
				sim = Similarity.D2(queryFeature, otherFeatures[i]);
			} else if (dis.equals("D2S")) {
				sim = Similarity.D2S(queryFeature, otherFeatures[i]);
			} else if (dis.equals("D2star")) {
				sim = Similarity.D2star(queryFeature, otherFeatures[i]);
			} else {
				sim = 0.0;
				System.out.println("dis输入有误");
			}
			genePair[i].setTotal_sim(sim);
			genePairs.add(genePair[i]);
			result.put(other[i], sim);
		}

		// 将map.entrySet()转换成list
		List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(result.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String, Double>>() {
			// 降序排序
			@Override
			public int compare(Entry<String, Double> o1, Entry<String, Double> o2) {
				// return o1.getValue().compareTo(o2.getValue());
				return o2.getValue().compareTo(o1.getValue()); // 从大到小排列，相似度越大，越相似
			}
		});

		int count1 = 0;
		String[] resultS = new String[39];

		for (Map.Entry<String, Double> mapping : list) {
			String tmp = count1 + ":" + mapping.getKey() + ":" + mapping.getValue();
			resultS[count1] = tmp;
			count1++;
			// System.out.println(tmp);
		}
		double[] acc = accuracy(resultS);
		double[] rs = new double[5];
		rs[0] = acc[0];
		rs[1] = acc[1];
		rs[2] = acc[2];
		rs[3] = acc[3];
		rs[4] = calROCRes(genePairs);
		if (code == 1 || code == 2 || code == 4) {
			double[][] coodinate = calROCCood(genePairs);
			try {
				WriteFile.writeCood(coodinate, "Sig" + dis + "-" + k + "-" + r);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return rs;
	}

	/**
	 * ROC计算
	 * 
	 * @param pos_pairs
	 * @param neg_pairs
	 * @return
	 */

	public static double calROCRes(List<GenePair> genepairs) {
		List<GenePair> negPairs = new ArrayList<>();
		List<GenePair> posPairs = new ArrayList<>();
		for (GenePair genePair : genepairs) {
			if (genePair.getName().indexOf("+") != -1) {
				posPairs.add(genePair);
			} else {
				negPairs.add(genePair);
			}
		}
		Collections.sort(genepairs, new Comparator<GenePair>() {
			public int compare(GenePair o1, GenePair o2) {
				if (o1.getTotal_sim() == o2.getTotal_sim())
					return 0;
				return o1.getTotal_sim() > o2.getTotal_sim() ? -1 : 1; // 应该从大到小
			}
		});

		int TP = 0, FP = 0, TN = negPairs.size(), FN = posPairs.size();
		double x = 0, y = 0;
		double area = 0;
		for (int i = 0; i < genepairs.size(); i++) {
			if (genepairs.get(i).isPos()) {
				TP++;
				FN--;
			} else {
				FP++;
				TN--;
			}
			area += (y + TP * 1.0 / (TP + FN)) * (FP * 1.0 / (FP + TN) - x) / 2.0;
			x = FP * 1.0 / (FP + TN);
			y = TP * 1.0 / (TP + FN);
		}
		// System.out.println("The AUC value is: " + area);
		return area;
	}

	/**
	 * 计算ROC坐标点
	 * 
	 * @param genepairs
	 * @param name
	 * @return
	 */
	public static double[][] calROCCood(List<GenePair> genepairs) {
		List<GenePair> negPairs = new ArrayList<>();
		List<GenePair> posPairs = new ArrayList<>();
		for (GenePair genePair : genepairs) {
			if (genePair.getName().indexOf("+") != -1) {
				posPairs.add(genePair);
			} else {
				negPairs.add(genePair);
			}
		}
		Collections.sort(genepairs, new Comparator<GenePair>() {
			public int compare(GenePair o1, GenePair o2) {
				if (o1.getTotal_sim() == o2.getTotal_sim())
					return 0;
				return o1.getTotal_sim() > o2.getTotal_sim() ? -1 : 1; // 应该从大到小
			}
		});

		int TP = 0, FP = 0, TN = negPairs.size(), FN = posPairs.size();
		double[][] cood = new double[genepairs.size()][2];
		for (int i = 0; i < genepairs.size(); i++) {
			if (genepairs.get(i).isPos()) {
				TP++;
				FN--;
			} else {
				FP++;
				TN--;
			}
			double TPR = TP * 1.0 / (TP + FN);
			double FPR = FP * 1.0 / (FP + TN);
			cood[i][0] = FPR;
			cood[i][1] = TPR;
		}
		return cood;
	}

	/**
	 * 熵值法加权
	 * 
	 * @param flag
	 * @param startk
	 * @param kmax
	 * @param code
	 * @param dis
	 * @param r
	 * @return
	 */
	public static double[] testAllkWeight2(boolean flag, int startk, int kmax, int code, String dis, int r) {
		// 读取目录
		ArrayList<String> filename = ReadFile.readAllFile("src/dataset/");
		String querySequence = ReadFile.readFileData("src/dataset/HSLIPAS.fasta");
		//
		querySequence = ReadFile.standData(querySequence);
		String[] other = new String[39];
		double[][] otherFeatures = new double[39][];
		double[][] otherfreqs = new double[39][];
		/**
		 * 相似度 结果result
		 */
		HashMap<String, Double> result = new HashMap<>();
		HashMap<String, String> map = new HashMap<>();
		for (int i = 0; i < filename.size(); i++) {
			String tmp1 = filename.get(i).split("\\.")[0].trim();
			if (tmp1.equals("HSLIPAS")) {

			} else {
				String tmp = "src/dataset/" + filename.get(i);
				String sequence = ReadFile.readFileData(tmp);
				map.put(filename.get(i).split("\\.")[0], sequence);
			}
		}
		// 获取查询序列的特征向量
		double[] queryFeature;
		// 查询序列所有频率
		double[] queryfreq = SeqComp.seqAllFeatureCount2(querySequence, startk, kmax, flag);
		if (code == 0) {
			queryFeature = SeqComp.seqAllFeaturefreq2(querySequence, startk, kmax, flag);
		} else if (code == 1) {
			queryFeature = SeqComp.seqAllFeatureCount2(querySequence, startk, kmax, flag);
		} else if (code == 2) {
			queryFeature = SeqComp.D2SAllCount2(querySequence, startk, kmax, r);
		} else if (code == 3) {
			queryFeature = SeqComp.haoAllFeature2(querySequence, startk, kmax);
		} else if (code == 4) {
			queryFeature = SeqComp.D2starAllCount2(querySequence, startk, kmax, r);
		} else {
			queryFeature = null;
			System.out.println("查询序列错误");
		}

		//
		/**
		 * 获取其余序列的特征向量
		 */
		Iterator<Map.Entry<String, String>> it = map.entrySet().iterator();
		int count = 0;
		String[] others = new String[39];

		ArrayList<GenePair> genePairs = new ArrayList<>();

		GenePair[] genePair = new GenePair[39];

		while (it.hasNext()) {
			Map.Entry<String, String> entry = it.next();
			boolean pos = false;
			if (entry.getKey().indexOf("+") != -1) {
				pos = true;
			} else {
				pos = false;
			}

			String seq = ReadFile.standData(entry.getValue());
			if (code == 0) {
				otherFeatures[count] = SeqComp.seqAllFeaturefreq2(seq, startk, kmax, flag);
			} else if (code == 1) {
				otherFeatures[count] = SeqComp.seqAllFeatureCount2(seq, startk, kmax, flag);
			} else if (code == 2) {
				otherFeatures[count] = SeqComp.D2SAllCount2(seq, startk, kmax, r);
			} else if (code == 3) {
				otherFeatures[count] = SeqComp.haoAllFeature2(seq, startk, kmax);
			} else if (code == 4) {
				otherFeatures[count] = SeqComp.D2starAllCount2(seq, startk, kmax, r);
			} else {
				otherFeatures[count] = null;
				System.out.println("查询序列错误");
			}
			otherfreqs[count] = SeqComp.seqAllFeaturefreq2(seq, startk, kmax, flag);
			others[count] = seq;
			other[count] = entry.getKey();
			int length = querySequence.length() + seq.length();
			// 构造基因对
			genePair[count] = new GenePair();
			genePair[count].setName(entry.getKey());
			genePair[count].setLength(length);
			genePair[count].setPos(pos);
			genePair[count].setSequence1(querySequence);
			genePair[count].setSequence2(seq);
			count++;
		}
		// 加上权重
		double[][] tmpFeature = new double[40][];
		tmpFeature[0] = queryFeature;
		double[][] freqs = new double[40][];
		freqs[0] = queryfreq;
		for (int i = 1; i < tmpFeature.length; i++) {
			tmpFeature[i] = otherFeatures[i - 1];
			freqs[i] = otherfreqs[i - 1];
		}
		// 使用freqs 计算权重
		double[] weight = SeqComp.getWeight1(freqs);
		double[][] pcaFeature = tmpFeature;

		// 计算相似度
		for (int i = 0; i < pcaFeature.length - 1; i++) {
			double sim = 0.0;
			if (dis.equals("eur")) {
				// sim = Similarity.euclidWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.euclid(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("cosine")) {
				sim = Similarity.cosine(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("ch")) {
				// sim = Similarity.chebyshevWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.chebyshev(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("ma")) {
				// sim = Similarity.manhattanWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.manhattan(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("pcc")) {
				sim = Similarity.pcc(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("kld")) {
				// sim = Similarity.KLDWeighted(pcaFeature[0], pcaFeature[i + 1],weight);
				sim = Similarity.KLD(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("hao")) {
				sim = Similarity.hao(pcaFeature[0], pcaFeature[i + 1]);
			} else if (dis.equals("D2")) {
				sim = Similarity.D2Weighted(pcaFeature[0], pcaFeature[i + 1], weight);
			} else if (dis.equals("D2S")) {
				sim = Similarity.D2SWeighted(pcaFeature[0], pcaFeature[i + 1], weight);
			} else if (dis.equals("D2star")) {
				sim = Similarity.D2starWeighted(queryFeature, otherFeatures[i], weight);
			} else {
				sim = 0.0;
				System.out.println("dis输入有误");
			}
			genePair[i].setTotal_sim(sim);
			genePairs.add(genePair[i]);
			result.put(other[i], sim);
		}

		// 将map.entrySet()转换成list
		List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(result.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String, Double>>() {
			// 降序排序
			@Override
			public int compare(Entry<String, Double> o1, Entry<String, Double> o2) {
				// return o1.getValue().compareTo(o2.getValue());
				return o2.getValue().compareTo(o1.getValue()); // 从大到小排列，相似度越大，越相似
			}
		});

		int count1 = 0;
		String[] resultS = new String[39];

		for (Map.Entry<String, Double> mapping : list) {
			String tmp = count1 + ":" + mapping.getKey() + ":" + mapping.getValue();
			resultS[count1] = tmp;
			count1++;
			// System.out.println(tmp);
		}
		double[] acc = accuracy(resultS);
		double[] rs = new double[5];
		rs[0] = acc[0];
		rs[1] = acc[1];
		rs[2] = acc[2];
		rs[3] = acc[3];
		rs[4] = calROCRes(genePairs);
		System.out.println("the  " + dis + " multik is:");
		return rs;
	}
}
