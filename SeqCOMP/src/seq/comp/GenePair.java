package seq.comp;
/**
 * »ùÒò¶Ô
 * @author zy
 *
 */
public class GenePair {
	private String name;
	private boolean pos;
	private String sequence1;
	private String sequence2;
	private int length;
	private double total_sim;
	
	public GenePair(String name,boolean pos, String sequence1, String sequence2,int length, double total_sim) {
		super();
		this.name=name;
		this.pos = pos;
		this.sequence1= sequence1;
		this.sequence2= sequence2;
		this.length = length;
		this.total_sim = total_sim;
	}
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public GenePair() {
		super();
	}

	public double getTotal_sim() {
		return total_sim;
	}
	public void setTotal_sim(double total_sim) {
		this.total_sim = total_sim;
	}
	public boolean isPos() {
		return pos;
	}
	public void setPos(boolean pos) {
		this.pos = pos;
	}
	public String getSequence1() {
		return sequence1;
	}
	public void setSequence1(String sequence) {
		this.sequence1 = sequence;
	}
	public String getSequence2() {
		return sequence2;
	}
	public void setSequence2(String sequence) {
		this.sequence2 = sequence;
	}
	public int getLength() {
		return length;
	}
	public void setLength(int length) {
		this.length = length;
	}
	

}
