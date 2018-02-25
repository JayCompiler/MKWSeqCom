package seq.comp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class ReadFile {
	private static ArrayList<String> listname = new ArrayList<String>();

	public static void main(String[] args) throws Exception {		
		String data1=readFileData("src/dataset/AGGGLINE-.fasta");
		System.out.println(data1.length());
		String data=standData(data1);
		StringBuffer sBuffer=new StringBuffer();
		for(int i=0;i<data.length();i++) {
			sBuffer.append(data.substring(i, i+1));
			if((i+1)%50==0) {
				System.out.println(sBuffer.toString());
				sBuffer.delete(0, sBuffer.length());
			}
		}
		System.out.println(data.length());
		System.out.println(data.indexOf("NNNNN"));

	}
	/**
	 * 读取一个目录下所有文件
	 * @param filepath
	 */
	public static ArrayList<String> readAllFile(String filepath) {
		File file = new File(filepath);
		String[] filelist = file.list();
		for (int i = 0; i < filelist.length; i++) {
			File readfile = new File(filepath + "/" + filelist[i]);
			listname.add(readfile.getName());
		}
		return listname;
	}
	/**
	 * 读取单个文件的数据
	 * @param filePath
	 * @param num
	 * @return
	 */
	public static String[]  readFileData2(String filePath,int num) {
		StringBuffer sBuffer=new StringBuffer();
		String tmp =new String();
		String [] data=new String[num];
		int count=0;
		try {
			BufferedReader br=new BufferedReader(new FileReader(filePath));
			try {
				while((tmp=br.readLine())!=null) {
					if(tmp.charAt(0)=='>'&&count==0) {
						count++;
						continue;
					}else if(tmp.charAt(0)=='>'&&count!=0) {	
						data[count-1]=standData(sBuffer.toString());
						sBuffer.delete(0,sBuffer.length());
						count++;
					}else {
						sBuffer.append(tmp);
					}
				}
				data[count-1]=standData(sBuffer.toString());
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return data;
	}
	
	/**
	 * 读取单个文件数据的数据名
	 * @param filePath
	 * @param num
	 * @return
	 */
	public static String[]  readFileDataName(String filePath,int num) {
		String tmp =new String();
		String[] name=new String[num];
		int count=0;
		try {
			BufferedReader br=new BufferedReader(new FileReader(filePath));
			try {
				while((tmp=br.readLine())!=null) {
					if(tmp.charAt(0)=='>'&&count==0) {
						name[count]=tmp.substring(1, tmp.length());
						count++;
					}else if(tmp.charAt(0)=='>'&&count!=0) {	
						name[count]=tmp.substring(1, tmp.length());
						count++;
					}else {
						continue;
					}
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return name;
	}
	/**
	 * 读取文件数据
	 * @param filePath
	 * @return
	 */
	public static String  readFileData(String filePath) {
		StringBuffer sBuffer=new StringBuffer();
		String tmp =new String();
		try {
			BufferedReader br=new BufferedReader(new FileReader(filePath));
			try {
				while((tmp=br.readLine())!=null) {
					//System.out.println((count++)+":"+tmp+"    "+tmp.charAt(0));
					if(tmp.charAt(0)!='>') 
						sBuffer.append(tmp);
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String data=sBuffer.toString();
		return data;
	}
	/**
	 * 去除字符串data中的非"A","T","C","G" 字符
	 * @param data
	 * @return
	 */
	public static String standData(String data) {
		int count=0;
		String tmpdata=data;
		for(int i=0;i<data.length();i++) {
			char tmp=data.charAt(i);
			if((tmp!='A')&&(tmp!='T')&&(tmp!='C')&&(tmp!='G')) {
				tmpdata=removeChar(i-count, tmpdata);
				count++;
			}
		}
		return tmpdata;
	}
	/**
	 * 删除字符串str中索引为index的字符
	 * @param index
	 * @param Str
	 * @return
	 */
	public static String removeChar(int index,String Str){  
	    Str=Str.substring(0,index)+Str.substring(index+1,Str.length());//substring的取值范围是:[,)  
	    return Str;  
	}  
}
