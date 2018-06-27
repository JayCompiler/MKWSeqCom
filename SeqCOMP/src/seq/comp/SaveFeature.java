package seq.comp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class SaveFeature {

	public static void main(String[] args) throws IOException {
		int k = 2;
		boolean flag = true;
		//保存特征数据的文件名：featureFileName,标签的文件名：labelFileName;
		String featureFileName="dataset1k2";
		String labelFileName="label1";
		String[] allk_mer = SeqComp.getK_mer(k);
		// 读取目录
		ArrayList<String> filename = ReadFile.readAllFile("src/dataset/");
		// 查询序列数据
		String querySequence = ReadFile.readFileData("src/dataset/HSLIPAS.fasta");
		// 获取查询序列的特征向量
		double[] queryFeature = SeqComp.seqFeaturefreq(querySequence, k, allk_mer, flag);
		querySequence = ReadFile.standData(querySequence);
		// 序列名称信息
		String[] other = new String[39];
		double[][] otherFeatures = new double[40][(int) Math.pow(4, k)];
		HashMap<String, String> map = new HashMap<>();
		for (int i = 0; i < filename.size(); i++) {
			String tmp1 = filename.get(i).split("\\.")[0].trim();
			// System.out.println(tmp1);
			if (tmp1.equals("HSLIPAS")) {

			} else {
				String tmp = "src/dataset/" + filename.get(i);
				String sequence = ReadFile.readFileData(tmp);
				map.put(filename.get(i).split("\\.")[0], sequence);
			}
		}
		/**
		 * 获取其余序列的特征向量
		 */
		Iterator<Map.Entry<String, String>> it = map.entrySet().iterator();
		int count = 0;
		while (it.hasNext()) {
			Map.Entry<String, String> entry = it.next();
			String pos = "";
			if (entry.getKey().indexOf("+") != -1) {
				pos = "1";
			} else {
				pos = "0";
			}
			String seq = ReadFile.standData(entry.getValue());
			otherFeatures[count] = SeqComp.seqFeaturefreq(seq, k, allk_mer, flag);
			other[count] = pos; //+ ":" + entry.getKey();
			count++;
		}
		//最后一条为查询序列特征
		otherFeatures[39]=queryFeature;
		for(int i=0;i<otherFeatures.length;i++) {
			otherFeatures[i]=SeqComp.normdata(otherFeatures[i]);
		}
//		double[][] positive =new double[20][];
//		
//		for(int i=0;i<other.length;i++) {
//			
//		}
		WriteFile.writeFeature(otherFeatures, featureFileName);
		WriteFile.writeLabel(other, labelFileName);
	}
}
