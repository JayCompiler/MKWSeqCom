package seq.comp;

import java.util.Random;

public class Experiment2 {
	/**
	 * test datasets
	 * "code" represents the pattern of extracting features:count ，freq， D2count haoFeature 
	 * "dis" represents the pattern of distance：eur,cosine,ch,ma
	 * pcc，kld,hao,D2，D2S，D2star
	 * 
	 * eur,cosine,ch,ma ，pcc (It doesn't have to be smooth) : flag=false code=0(freq); 
	 * kld, smoothening required flag=true code=0
	 * D2 flag =flase ,code=1 (count) 
	 * D2S，smoothening required flag=true code=2(D2count)  r=0 or 1 or 2
	 * hao smoothening required flag=true code =3 
	 * D2star smoothening required  flag=true code=4(D2starcount)  r=0 or 1 or 2
	 */
	public static void test() {
		int k=7;
		int r=1;
		String dis="pcc";
		int code=0;
		boolean flag=true;
		// 读取数据
		//String[] data = ReadFile.readFileData2("src/fRS/fly_blastoderm.fa", 164);
		//String[] data = ReadFile.readFileData2("src/fRS/fly_eye.fa", 34);
		//String[] data = ReadFile.readFileData2("src/fRS/fly_pns.fa", 46);
		//String[] data = ReadFile.readFileData2("src/fRS/fly_tracheal_system.fa", 18);
		//String[] data = ReadFile.readFileData2("src/fRS/human_HBB.fa", 34);
		//String[] data = ReadFile.readFileData2("src/fRS/human_liver.fa", 18);
		String[] data = ReadFile.readFileData2("src/fRS/human_muscle.fa", 56);
		int size = data.length;
		int halfsize=size/2;
		String[] pos = new String[halfsize];
		String[] neg = new String[halfsize];
		//赋值初始化
		for(int i=0;i<halfsize;i++) {
			pos[i]=data[i];
			neg[i]=data[i+halfsize];
		}
		int pairsize = halfsize * (halfsize - 1) / 2;
		GenePair[] posPairs = new GenePair[pairsize];
		GenePair[] negPairs = new GenePair[pairsize];
		int count=0;
		//生成基因对
		for (int i = 1; i < halfsize; i++) {
			for (int j = i + 1; j <= halfsize; j++) {
				int length1=pos[i-1].length()+pos[j-1].length();
				int length2=neg[i-1].length()+neg[j-1].length();
				posPairs[count]=new GenePair(i+"-"+j, true, pos[i-1], pos[j-1], length1, 0.0);
				negPairs[count]=new GenePair(-i+"-"+j, false, neg[i-1], neg[j-1], length2, 0.0);
				count++;
			}
		}
		/*计算相似度*/
		for (int i=0;i<posPairs.length;i++) {
			posPairs[i].setTotal_sim(computeSim(posPairs[i], dis, k, code, r, flag));
			negPairs[i].setTotal_sim(computeSim(negPairs[i], dis, k, code, r, flag));
		}
		GenePair[] all=mergeSortGenePair(posPairs, negPairs);
		int top=(halfsize-1)*halfsize/2;
		if(top>300) {
			all=select(all, 300);
		}else {
			all=select(all, top);
		}

		double [] ratio =getRatio(all);
		for (double d : ratio) {
			System.out.print(d+" ");
		}
		System.out.println();
		int [] num=getNum(all);
		for (int i : num) {
			System.out.print(i+" ");
		}
	}
	
	/**
	 * test multiply k-value weighted
	 * D2 flag =flase ,code=1 (count) 
	 * D2S， smoothening required flag=true code=2(D2count)  r=0 or 1 or 2
	 * D2star smoothening required flag=true code=4(D2starcount)  r=0 or 1 or 2
	 */
	public static void testMultiK() {
		int kmax=7;
		String dis="D2star";
		int code=4;
		int r=1;
		int startK=2;
		boolean flag=true;
		// 读取数据
		//String[] data = ReadFile.readFileData2("src/fRS/fly_blastoderm.fa", 164);
		//String[] data = ReadFile.readFileData2("src/fRS/fly_eye.fa", 34);
		//String[] data = ReadFile.readFileData2("src/fRS/fly_pns.fa", 46);
		String[] data = ReadFile.readFileData2("src/fRS/fly_tracheal_system.fa", 18);
		//String[] data = ReadFile.readFileData2("src/fRS/human_HBB.fa", 34);
		//String[] data = ReadFile.readFileData2("src/fRS/human_liver.fa", 18);
		//String[] data = ReadFile.readFileData2("src/fRS/human_muscle.fa", 56);
		int size = data.length;
		int halfsize=size/2;
		String[] pos = new String[halfsize];
		String[] neg = new String[halfsize];
		//赋init
		for(int i=0;i<halfsize;i++) {
			pos[i]=data[i];
			neg[i]=data[i+halfsize];
		}
		int pairsize = halfsize * (halfsize - 1) / 2;
		GenePair[] posPairs = new GenePair[pairsize];
		GenePair[] negPairs = new GenePair[pairsize];
		double[][] freqs=new double[size][];
		
		for(int i=0;i<size;i++) {
			freqs[i]=SeqComp.seqAllFeaturefreq(data[i],startK, kmax, true);
		}
		//computer weight
		double [] Weight=SeqComp.getWeight(freqs);
		int count=0;
		//generate gene_pairs
		for (int i = 1; i < halfsize; i++) {
			for (int j = i + 1; j <= halfsize; j++) {
				int length1=pos[i-1].length()+pos[j-1].length();
				int length2=neg[i-1].length()+neg[j-1].length();
				posPairs[count]=new GenePair(i+"-"+j, true, pos[i-1], pos[j-1], length1, 0.0);
				negPairs[count]=new GenePair(-i+"-"+j, false, neg[i-1], neg[j-1], length2, 0.0);
				count++;
			}
		}
		/*computer similarity*/
		for (int i=0;i<posPairs.length;i++) {
			posPairs[i].setTotal_sim(computeMultKSim(posPairs[i], dis,startK, kmax, code, r, flag, Weight));
			negPairs[i].setTotal_sim(computeMultKSim(negPairs[i], dis, startK,kmax, code, r, flag, Weight));
		}
		GenePair[] all=mergeSortGenePair(posPairs, negPairs);
		int top=halfsize*(halfsize-1)/2;
		if(top>300) {
			all=select(all, 300);
		}else {
			all=select(all, top);
		}
		double [] ratio =getRatio(all);
		for (double d : ratio) {
			System.out.print(d+" ");
		}
		System.out.println();
		int [] num=getNum(all);
		for (int i : num) {
			System.out.print(i+" ");
		}
	}
	/**
	 * Merge String[]
	 * @param s1
	 * @param s2
	 * @return
	 */
	public static String [] merSt(String[] s1,String[] s2) {
		int len1=s1.length;
		int len2=s2.length;
		int len=len1+len2;
		String[] string=new String[len];
		int count=0;
		for (int i=0;i<len1;i++) {
			string[count++]=s1[i];
		}
		for (int i=0;i<len2;i++) {
			string[count++]=s2[i];
		}
		return string;
	}
	/**
	 * Random choose n sequeces from seqs
	 * @param seqs
	 * @param n
	 * @return
	 */
	public static String[] shuffle(String[] seqs,int n) {
		String [] select=new String[n];
		int [] num=randomArray(0, seqs.length-1,n);
		for(int i=0;i<num.length;i++) {
			select[i]=seqs[num[i]];
		}
		return select;
	}
	/**
	 * Generate n different integers in a range
	 * @param min
	 * @param max
	 * @param n
	 * @return
	 */
	public static int[] randomArray(int min,int max,int n){  
	    int len = max-min+1;  
	      
	    if(max < min || n > len){  
	        return null;  
	    }  
	      
	    //init
	    int[] source = new int[len];  
	       for (int i = min; i < min+len; i++){  
	        source[i-min] = i;  
	       }  
	         
	       int[] result = new int[n];  
	       Random rd = new Random();  
	       int index = 0;  
	       for (int i = 0; i < result.length; i++) {  
	           index = Math.abs(rd.nextInt() % len--);  
	           result[i] = source[index];  
	           source[index] = source[len];  
	       }  
	       return result;  
	}  

	/**
	 * compute rate
	 * @param genePairs
	 * @return
	 */
	public static double[] getRatio(GenePair[] genePairs) {
		int size=genePairs.length;
		int m=(int)Math.ceil((double)size/10.0);
		double [] ratio=new double[m];
		int count=0;
		for(int i=0;i<m;i++) {
			int all=0;
			if(i==(m-1)) {
				for(int j=10*i;j<size;j++) {
					if(genePairs[j].getName().charAt(0)!='-') {
						count++;
					}
				}
				all=size;
			}else {
				for(int j=10*i;j<10*(i+1);j++) {
					if(genePairs[j].getName().charAt(0)!='-') {
						count++;
					}
				}
				all=10*(i+1);
			}
			ratio[i]=(double)count/(double)all;
		}
		return ratio;
	}
	/**
	 * computer number
	 * @param genePairs
	 * @return
	 */
	public static int[] getNum(GenePair[] genePairs) {
		int size=genePairs.length;
		int m=(int)Math.ceil((double)size/10.0);
		int [] num=new int[m];
		int count=0;
		for(int i=0;i<m;i++) {
			if(i==(m-1)) {
				for(int j=10*i;j<size;j++) {
					if(genePairs[j].getName().charAt(0)!='-') {
						count++;
					}
				}
			}else {
				for(int j=10*i;j<10*(i+1);j++) {
					if(genePairs[j].getName().charAt(0)!='-') {
						count++;
					}
				}
			}
			num[i]=count;
		}
		return num;
	}
	/**
	 * computer multiple k-value: D2,D2S,D2star
	 * @param X
	 * @param dis
	 * @param k
	 * @param code
	 * @param r
	 * @param flag
	 * @return
	 */
	public static double computeMultKSim(GenePair X,String dis,int startK,int kmax,int code,int r,boolean flag,double[] weight1) {
		int size=0;
		//提取特征
		for(int i=2;i<=kmax;i++) {
			size=(int)Math.pow(4, i)+size;
		}
		double[] feature1=new double[size];
		double[] feature2=new double[size];
		double[] freq1=SeqComp.seqAllFeaturefreq(X.getSequence1(),startK, kmax, flag);
		double[] freq2=SeqComp.seqAllFeaturefreq(X.getSequence2(),startK, kmax, flag);
		double [][] freqs=new double[2][];
		freqs[0]=freq1;
		freqs[1]=freq2;
		if(code==1) {
			feature1 = SeqComp.seqAllFeatureCount(X.getSequence1(),startK, kmax, flag);
			feature2 = SeqComp.seqAllFeatureCount(X.getSequence2(),startK, kmax, flag);
		}else if(code==2) {
			feature1 = SeqComp.D2SAllCount(X.getSequence1(),startK, kmax, r);
			feature2 = SeqComp.D2SAllCount(X.getSequence2(), startK,kmax, r);
		}else if(code==4) {
			feature1 = SeqComp.D2starAllCount(X.getSequence1(),startK, kmax, r);
			feature2 = SeqComp.D2starAllCount(X.getSequence2(),startK, kmax, r);
		}else {
			System.out.println("input code is error,only support 1，2，4");
		}
		//计算相似度
		double sim=0.0;
		if (dis.equals("D2")) {
			sim = Similarity.D2Weighted(feature1, feature2,weight1);
		}  else if (dis.equals("D2S")) {
			sim = Similarity.D2SWeighted(feature1, feature2,weight1);
		} else if (dis.equals("D2star")) {
			sim = Similarity.D2starWeighted(feature1, feature2,weight1);
		} else {
			System.out.println("the 'dis' is error");
		}
		return sim;
	}
	/**
	 * computer the similarity of gene sequences
	 * */
	public static double computeSim(GenePair X,String dis,int k,int code,int r,boolean flag) {
		String[] allk_mer=SeqComp.getK_mer(k);
		//提取特征
		int size=(int)Math.pow(4, k);
		double[] feature1=new double[size];
		double[] feature2=new double[size];
		if(code==0) {
			feature1 = SeqComp.seqFeaturefreq(X.getSequence1(), k, allk_mer, flag);
			feature2 = SeqComp.seqFeaturefreq(X.getSequence2(), k, allk_mer, flag);
		}else if(code==1) {
			feature1 = SeqComp.seqFeatureCount(X.getSequence1(), k, allk_mer, flag);
			feature2 = SeqComp.seqFeatureCount(X.getSequence2(), k, allk_mer, flag);
		}else if(code==2) {
			feature1 = SeqComp.D2SCount(X.getSequence1(), k, allk_mer, r);
			feature2 = SeqComp.D2SCount(X.getSequence2(), k, allk_mer, r);
		}else if(code==3) {
			feature1 = SeqComp.haoFeature(X.getSequence1(), k);
			feature2 = SeqComp.haoFeature(X.getSequence2(), k);
		}else if(code==4) {
			feature1 = SeqComp.D2starCount(X.getSequence1(), k, allk_mer, r);
			feature2 = SeqComp.D2starCount(X.getSequence2(), k, allk_mer, r);
		}
		double sim=0.0;
		if (dis.equals("eur")) {
			sim = Similarity.euclid(feature1, feature2);
		} else if (dis.equals("cosine")) {
			sim = Similarity.cosine(feature1, feature2);
		} else if (dis.equals("ch")) {
			sim = Similarity.chebyshev(feature1,feature2);
		} else if (dis.equals("ma")) {
			sim = Similarity.manhattan(feature1, feature2);
		} else if (dis.equals("pcc")) {
			sim = Similarity.pcc(feature1,feature2);
		} else if (dis.equals("kld")) {
			sim = Similarity.KLD(feature1, feature2);
		} else if (dis.equals("hao")) {
			sim = Similarity.hao(feature1, feature2);
		} else if (dis.equals("D2")) {
			sim = Similarity.D2(feature1, feature2);
		} else if (dis.equals("D2S")) {
			sim = Similarity.D2S(feature1, feature2);
		} else if (dis.equals("D2star")) {
			sim = Similarity.D2star(feature1, feature2);
		} else {
			System.out.println("the 'dis' is error");
		}
		return sim;
	}
	
	/**
	 * According to the order of similarity,
		sort and merge, the bigger and the more similar
	 * @param poss
	 * @param negs
	 * @return
	 */
	public static GenePair[] mergeSortGenePair(GenePair[] poss,GenePair[] negs) {
		int length=poss.length+negs.length;
		GenePair[]  genePairs=new GenePair[length];
		int count=0;
		for(int i=0;i<poss.length;i++) {
			genePairs[count++]=poss[i];
		}
		for(int i=0;i<poss.length;i++) {
			genePairs[count++]=negs[i];
		}
		sortGenePairs(genePairs);
		return genePairs;
	}
	/**
	 * sort
	 * @param genepairs
	 */
	public static void sortGenePairs(GenePair[] genepairs) {
		int size=genepairs.length;
		for(int i=0;i<size-1;i++) {
			for(int j=i+1;j<size;j++) {
				if(genepairs[i].getTotal_sim()<genepairs[j].getTotal_sim()) {
					swap(genepairs, i, j);
				}
			}
		}
	}
	/**
	 * change
	 * @param genepairs
	 * @param i
	 * @param j
	 */
	private static void swap(GenePair[] genepairs,int i,int j) {
		GenePair tmPair=genepairs[i];
		genepairs[i]=genepairs[j];
		genepairs[j]=tmPair;
	}
	/**
	 * Select the first k gene pairs in the genePairs sequence
	 * @param genePairs
	 * @param k
	 * @return
	 */
	public static GenePair[] select(GenePair[] genePairs,int k) {
		GenePair[] data=new GenePair[k];
		for(int i=0;i<k;i++) {
			data[i]=genePairs[i];
		}
		return data;
	}
}
