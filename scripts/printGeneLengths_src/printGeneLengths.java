package printGeneLengths;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;

public class printGeneLengths {

//	Required arguments:
//	1. Annotation file in GTF format
//	2. Path to output file
	public static void main(String[] args) throws IOException {
//		Read in exons of every gene from GTF annotation file
		HashMap<String, GTFEntry> gtffile = GTFHelper.readGTF(args[0]);

//		Write out length of exonic region for every gene
		BufferedWriter bw = new BufferedWriter(new FileWriter(args[1]));
		bw.append("chr\tstart\tend\tgeneID\tgeneName\tlength\n");
		for(Entry<String,GTFEntry> e : gtffile.entrySet()) {
			bw.append(e.getValue().getChr()+"\t"+e.getValue().getStart()+"\t"+e.getValue().getEnd()+"\t"+e.getKey()+"\t"+e.getValue().getGeneName()+"\t"+e.getValue().getGeneLength()+"\n");
			bw.flush();
		}
		bw.close();
	}

}
