package printGeneLengths;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class GTFHelper {

//	Parses exonic regions for every gene in GTF file
//	Overlapping exonic regions are counted only once
	public static HashMap<String, GTFEntry> readGTF(String file) {
		HashMap<String, GTFEntry> result = new HashMap<String, GTFEntry>();
		String line;
		String[] temp;
		String[] temp2;
		String entry;
		String genename;
		String chr;
		int start;
		int end;

		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
				if(line.startsWith("#")) continue;
				temp = line.split("\t");

				if(!temp[2].equals("exon") && !temp[2].equals("gene")) continue;

				chr = temp[0];
				start = Integer.parseInt(temp[3]);
				end = Integer.parseInt(temp[4]);
				temp2 = temp[8].split("gene_id \"");
				temp2 = temp2[1].split("\"");
				entry = temp2[0];
				
				if(temp[2].equals("gene")) {
					temp2 = temp[8].split("gene_name \"");
					temp2 = temp2[1].split("\"");
					genename = temp2[0];
					result.put(entry, new GTFEntry(chr,entry,genename,start,end));
				} else if(temp[2].equals("exon")){
					result.get(entry).addExon(start, end);
				} else {
					continue;
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return result;
	}
}
