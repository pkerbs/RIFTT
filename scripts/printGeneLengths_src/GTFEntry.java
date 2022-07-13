package printGeneLengths;

import java.util.LinkedList;

public class GTFEntry{

	private String chr;
	private int start;
	private int end;
	private String entry;
	private String genename;
	private LinkedList<Exon> exons;
	
	public GTFEntry(String chr, String entry, String genename,int start,int end){
		this.chr = chr;
		this.entry=entry;
		this.genename=genename;
		this.start=start;
		this.end=end;
		this.exons = new LinkedList<Exon>();
	}
	
	public String getGeneID(){
		return this.entry;
	}
	
	public String getGeneName(){
		return this.genename;
	}
	
	public String getChr(){
		return this.chr;
	}
	
	public int getStart(){
		return this.start;
	}
	
	public int getEnd(){
		return this.end;
	}
	
	public void addExon(int start, int end){
		int index=this.getIndexOfOverlappedExon(start, end);
		if(index > -1) {
			this.exons.get(index).setStart(Math.min(start, this.exons.get(index).getStart()));
			this.exons.get(index).setEnd(Math.max(end, this.exons.get(index).getEnd()));
		} else {
			this.exons.add(new Exon(start,end));
		}
	}
	
	public int getGeneLength(){
		int length = 0;
		for(Exon e:this.exons) {
			length += (e.getEnd()-e.getStart()+1);
		}
		return length;
	}
		
	public int getIndexOfOverlappedExon(int start, int end){
		for (int i = 0; i < this.exons.size(); i++) {
			if(this.exons.get(i).overlapsWith(start,end)) return i;
		}
		return -1;
	}
	
	public LinkedList<Exon> getExons(){
		return this.exons;
	}
	
	public class Exon{
		private int start;
		private int end;	
		
		public Exon(int start, int end){
			this.start=start;
			this.end=end;
		}
		
		public int getStart(){
			return this.start;
		}
		
		public int getEnd(){
			return this.end;
		}
		
		public void setStart(int start){
			this.start = start;
		}
		
		public void setEnd(int end){
			this.end = end;
		}
		
		public boolean overlapsWith(int start, int end){
			if((this.start > start && this.start > end) || (this.end < start && this.end < end)) return false;
			else return true;
		}
	}
}
