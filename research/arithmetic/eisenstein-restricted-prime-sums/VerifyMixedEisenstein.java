import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HexFormat;
import java.util.List;

/**
 * Independent exact verifier for N=a+3b, a,b supported on primes1mod3.
 * The integer1 is admitted. This program reads no Python data.
 *
 * java -Xmx512m VerifyMixedEisenstein.java 2000001 mixed_java_receipt.json
 *
 * The first argument is an exclusive bound. Indexj denotes6j+1 and
 * N=6k+4 is represented exactly when allowed[j] and allowed[k-3j].
 */
public class VerifyMixedEisenstein {
    public static void main(String[] args) throws Exception {
        long started=System.nanoTime();
        int limit=args.length>0?Integer.parseInt(args[0]):2000001;
        Path output=Path.of(args.length>1?args[1]:"mixed_java_receipt.json");
        if(limit<=4)throw new IllegalArgumentException("Exclusive bound must exceed4");
        int lastIndex=(limit-5)/6;
        int lastTarget=6*lastIndex+4;
        int maximumSummand=lastTarget-3;
        int root=(int)Math.sqrt(maximumSummand);
        while((long)(root+1)*(root+1)<=maximumSummand)root++;
        boolean[] composite=new boolean[root+1];
        for(int p=2;(long)p*p<=root;p++)if(!composite[p])
            for(int m=p*p;m<=root;m+=p)composite[m]=true;
        byte[] allowed=new byte[lastIndex+1];
        Arrays.fill(allowed,(byte)1);
        for(int p=5;p<=root;p+=6)if(!composite[p])
            for(long m=(long)p*p;m<=maximumSummand;m+=6L*p)
                allowed[(int)((m-1)/6)]=0;

        long[] histogram=new long[128];
        List<Integer> uncovered=new ArrayList<>();
        List<Integer> maximumTargets=new ArrayList<>();
        int maximumWitness=-1;
        MessageDigest witnessDigest=MessageDigest.getInstance("SHA-256");
        byte[] digestBuffer=new byte[65536];
        int bufferPosition=0;
        for(int k=0;k<=lastIndex;k++) {
            int witness=-1;
            for(int j=0;j<=k/3;j++) {
                if(allowed[j]==0||allowed[k-3*j]==0)continue;
                witness=j;
                if(j>=histogram.length)histogram=Arrays.copyOf(histogram,Math.max(j+1,2*histogram.length));
                histogram[j]++;
                if(j>maximumWitness){maximumWitness=j;maximumTargets.clear();}
                if(j==maximumWitness)maximumTargets.add(6*k+4);
                break;
            }
            if(witness<0)uncovered.add(6*k+4);
            digestBuffer[bufferPosition++]=(byte)witness;
            digestBuffer[bufferPosition++]=(byte)(witness>>>8);
            digestBuffer[bufferPosition++]=(byte)(witness>>>16);
            digestBuffer[bufferPosition++]=(byte)(witness>>>24);
            if(bufferPosition==digestBuffer.length){witnessDigest.update(digestBuffer);bufferPosition=0;}
        }
        witnessDigest.update(digestBuffer,0,bufferPosition);
        StringBuilder out=new StringBuilder("{\n");
        out.append("  \"scope\": \"Independent Java finite verification ofN=a+3b;1 allowed;no input witness files\",\n");
        out.append("  \"exclusive_limit\": ").append(limit).append(",\n");
        out.append("  \"largest_target\": ").append(lastTarget).append(",\n");
        out.append("  \"number_of_targets\": ").append(lastIndex+1).append(",\n");
        out.append("  \"exceptions\": ").append(uncovered).append(",\n");
        out.append("  \"uncovered_count\": ").append(uncovered.size()).append(",\n");
        out.append("  \"largest_least_b\": ").append(6*maximumWitness+1).append(",\n");
        out.append("  \"largest_least_weighted_summand_3b\": ").append(18*maximumWitness+3).append(",\n");
        out.append("  \"largest_least_b_targets\": ").append(maximumTargets).append(",\n");
        out.append("  \"least_b_histogram\": [\n");
        boolean first=true;
        for(int j=0;j<histogram.length;j++)if(histogram[j]!=0) {
            if(!first)out.append(",\n");
            first=false;
            out.append("    {\"b\": ").append(6*j+1).append(", \"weighted_summand_3b\": ").append(18*j+3)
                .append(", \"new_targets\": ").append(histogram[j]).append("}");
        }
        out.append("\n  ],\n");
        out.append("  \"sha256_admissibility_bits\": \"").append(HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(allowed))).append("\",\n");
        out.append("  \"sha256_int32_witness_indices\": \"").append(HexFormat.of().formatHex(witnessDigest.digest())).append("\",\n");
        out.append("  \"java_runtime\": \"").append(Runtime.version()).append("\",\n");
        Path source=Path.of("VerifyMixedEisenstein.java");
        if(Files.isRegularFile(source))out.append("  \"source_sha256\": \"").append(HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(source)))).append("\",\n");
        out.append("  \"elapsed_seconds\": ").append((System.nanoTime()-started)/1_000_000_000.0).append("\n}\n");
        Files.writeString(output,out);
        System.out.print(out);
        if(!uncovered.isEmpty())throw new IllegalStateException("Exact counterexamples found");
    }
}
