package seq.comp;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

public class Main1 {

	/**
	 * code������ȡ������ģʽ����count ��freq�� D2count haoFeature dis ��ʾ����ģʽ��eur,cosine,ch,ma
	 * ��pcc��kld,hao,D2��D2S��D2star
	 * 
	 * eur,cosine,ch,ma ��pcc ����Ҫƽ�� : flag=false code=0(freq); 
	 * kld, ��Ҫƽ�� flag=true code=0
	 * D2 flag =flase ,code=1 (count) 
	 * D2S����Ҫƽ�� flag=true code=2(D2count)  r=0��1
	 * hao���� ��Ҫƽ�� flag=true code =3 
	 * D2star��Ҫƽ�� flag=true code=4(D2starcount)  r=0��1
	 *  
	 *ȫ k pca ���� flag=flase code��ʱδ����eur����
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		HashMap<String, Integer> list = new HashMap<>();
		
		list.put("eur", 0);
		list.put("cosine", 0);
		list.put("ch", 0);
		list.put("ma", 0);
		list.put("pcc", 0);
		list.put("kld", 0);
		list.put("D2", 1);
		list.put("D2S", 2);
		list.put("D2star", 4);
		list.put("hao", 3);
		// int code=0;
		// String dis="eur";
		boolean flag = true; // �Ƿ���Ҫƽ��
		int r =0;
		 int kmax=6;
		 int startK=2;
		/* ������pca�������������kֵ */
//		 double[][] acc=new double[list.size()][5];
//		 Iterator<Entry<String, Integer>> iterator=list.entrySet().iterator();
//		 int count=0;
//		 while(iterator.hasNext()) {
//		 Entry<String, Integer> entry=iterator.next();
//		 acc[count]=Experiment1.testAllkWeight(flag,startK,kmax, entry.getValue(), entry.getKey(),r);
//		 //acc[count]=Experiment1.testAllkNoPCA(flag, startK, kmax, entry.getValue(), entry.getKey(), r);
//		 for (double ds : acc[count]) {
//		 System.out.print(ds+" ");
//		 }
//		 
//		 count++;
//		 System.out.println();
//		 }



		/* ����kֵ */
		double[][][] acc = new double[list.size()][7][5];
		Iterator<Entry<String, Integer>> iterator = list.entrySet().iterator();
		int count = 0;
		while (iterator.hasNext()) {
			Entry<String, Integer> entry = iterator.next();
			int count1 = 0;
			for (int k = 2; k <= 8; k++) {
				acc[count][count1] = Experiment1.test(k, r, flag, entry.getValue(), entry.getKey());
				count1++;
			}
			System.out.println("the " + entry.getKey() + " result is:");
			for (double[] ds : acc[count]) {
				System.out.println();
				for (double d : ds) {
					System.out.print(d + " ");
				}
				System.out.println();
			}
			count++;
			System.out.println();
		}

	}

}
