import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HexFormat;
import java.util.List;

/**
 * Independent exact finite verifier; requires Java 21 and no dependencies.
 *
 * From this directory:
 *   java VerifyEisenstein.java 100000000 java_finite_receipt.json
 *
 * The first argument is an EXCLUSIVE target bound. This program reads no
 * Python output or witness file. It constructs its own sieve and checks each
 * target directly against every allowed smaller summand up to 211.
 *
 * Index j represents 6j+1. Such integers already exclude 2 and 3. An excluded
 * prime q=5 (mod6) must occur with even total multiplicity across the excluded
 * prime factors. Consequently at least one excluded prime is <=sqrt(n).
 * Sieving from q*q in steps 6*q therefore removes precisely the forbidden
 * integers; smaller forbidden multiples were removed at a smaller prime.
 *
 * The balanced search independently examines candidates starting at N/2 and
 * descending to N/3. Neither search proves anything above the stated bound.
 */
public class VerifyEisenstein {
    private static final int MAX_SMALL_SUMMAND = 211;

    private static final class IntDigest {
        private final MessageDigest digest;
        private final byte[] buffer = new byte[65536];
        private int position;

        IntDigest() throws Exception {
            digest = MessageDigest.getInstance("SHA-256");
        }

        void append(int value) {
            // Match the documented little-endian signed int32 witness encoding.
            buffer[position++] = (byte) value;
            buffer[position++] = (byte) (value >>> 8);
            buffer[position++] = (byte) (value >>> 16);
            buffer[position++] = (byte) (value >>> 24);
            if (position == buffer.length) {
                digest.update(buffer);
                position = 0;
            }
        }

        String finish() {
            digest.update(buffer, 0, position);
            return HexFormat.of().formatHex(digest.digest());
        }
    }

    private static String hash(byte[] data) throws Exception {
        return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(data));
    }

    private static String jsonList(List<Integer> values) {
        return values.toString();
    }

    public static void main(String[] args) throws Exception {
        long started = System.nanoTime();
        int exclusiveLimit = args.length > 0 ? Integer.parseInt(args[0]) : 100000000;
        Path output = Path.of(args.length > 1 ? args[1] : "java_finite_receipt.json");
        if (exclusiveLimit <= 2) {
            throw new IllegalArgumentException("Exclusive limit must exceed 2");
        }

        int lastIndex = (exclusiveLimit - 3) / 6;
        int lastTarget = 6 * lastIndex + 2;
        int maximumSummand = lastTarget - 1;
        int sieveLimit = (int) Math.sqrt(maximumSummand);
        while ((long) (sieveLimit + 1) * (sieveLimit + 1) <= maximumSummand) sieveLimit++;
        while ((long) sieveLimit * sieveLimit > maximumSummand) sieveLimit--;

        boolean[] composite = new boolean[sieveLimit + 1];
        for (int p = 2; (long) p * p <= sieveLimit; p++) {
            if (!composite[p]) {
                for (int multiple = p * p; multiple <= sieveLimit; multiple += p) {
                    composite[multiple] = true;
                }
            }
        }

        byte[] allowed = new byte[lastIndex + 1];
        Arrays.fill(allowed, (byte) 1);
        int excludedPrimes = 0;
        for (int p = 5; p <= sieveLimit; p += 6) {
            if (composite[p]) continue;
            excludedPrimes++;
            for (long value = (long) p * p; value <= maximumSummand; value += 6L * p) {
                allowed[(int) ((value - 1) / 6)] = 0;
            }
        }
        int allowedCount = 0;
        for (byte member : allowed) allowedCount += member;

        List<Integer> candidates = new ArrayList<>();
        for (int a = 1; a <= MAX_SMALL_SUMMAND && a <= maximumSummand; a += 6) {
            if (allowed[(a - 1) / 6] != 0) candidates.add(a);
        }
        long[] histogram = new long[candidates.size()];
        List<Integer> uncovered = new ArrayList<>();
        List<Integer> balancedExceptions = new ArrayList<>();
        List<Integer> maximumLeastTargets = new ArrayList<>();
        IntDigest unrestrictedDigest = new IntDigest();
        IntDigest balancedDigest = new IntDigest();
        int maximumLeast = 0;
        int maximumBalancedOffset = 0;

        for (int k = 0; k <= lastIndex; k++) {
            int n = 6 * k + 2;
            int witness = -1;
            for (int c = 0; c < candidates.size(); c++) {
                int a = candidates.get(c);
                if (2L * a > n) break;
                int b = n - a;
                if (allowed[(b - 1) / 6] == 0) continue;
                witness = (a - 1) / 6;
                histogram[c]++;
                if (a > maximumLeast) {
                    maximumLeast = a;
                    maximumLeastTargets.clear();
                }
                if (a == maximumLeast) maximumLeastTargets.add(n);
                break;
            }
            unrestrictedDigest.append(witness);
            if (witness < 0) uncovered.add(n);

            int balancedWitness = -1;
            for (int t = k / 2; 3L * t >= k; t--) {
                if (allowed[t] != 0 && allowed[k - t] != 0) {
                    balancedWitness = t;
                    maximumBalancedOffset = Math.max(maximumBalancedOffset, k / 2 - t);
                    break;
                }
            }
            balancedDigest.append(balancedWitness);
            if (balancedWitness < 0) balancedExceptions.add(n);
        }

        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"scope\": \"Independent Java exact finite verification; no input witnesses\",\n");
        json.append("  \"exclusive_limit\": ").append(exclusiveLimit).append(",\n");
        json.append("  \"largest_target\": ").append(lastTarget).append(",\n");
        json.append("  \"number_of_targets\": ").append(lastIndex + 1).append(",\n");
        json.append("  \"allowed_summands_count\": ").append(allowedCount).append(",\n");
        json.append("  \"sieve_prime_limit\": ").append(sieveLimit).append(",\n");
        json.append("  \"odd_offending_primes_sieved\": ").append(excludedPrimes).append(",\n");
        json.append("  \"tested_small_summand_limit\": ").append(MAX_SMALL_SUMMAND).append(",\n");
        json.append("  \"uncovered_targets\": ").append(jsonList(uncovered)).append(",\n");
        json.append("  \"largest_least_summand\": ").append(maximumLeast).append(",\n");
        json.append("  \"largest_least_summand_targets\": ").append(jsonList(maximumLeastTargets)).append(",\n");
        json.append("  \"least_summand_histogram\": [\n");
        for (int c = 0; c < candidates.size(); c++) {
            json.append("    {\"summand\": ").append(candidates.get(c))
                .append(", \"new_targets\": ").append(histogram[c]).append("}")
                .append(c + 1 == candidates.size() ? "\n" : ",\n");
        }
        json.append("  ],\n");
        json.append("  \"balanced_condition\": \"N/3 <= a <= b <= 2N/3\",\n");
        json.append("  \"balanced_exceptions\": ").append(jsonList(balancedExceptions)).append(",\n");
        json.append("  \"balanced_max_search_offset_in_index_units\": ").append(maximumBalancedOffset).append(",\n");
        json.append("  \"sha256_allowed_boolean_array\": \"").append(hash(allowed)).append("\",\n");
        json.append("  \"sha256_unrestricted_int32_witness_indices\": \"").append(unrestrictedDigest.finish()).append("\",\n");
        json.append("  \"sha256_balanced_int32_witness_indices\": \"").append(balancedDigest.finish()).append("\",\n");
        json.append("  \"java_runtime\": \"").append(Runtime.version()).append("\",\n");
        Path source = Path.of("VerifyEisenstein.java");
        if (Files.isRegularFile(source)) {
            json.append("  \"source_sha256\": \"").append(hash(Files.readAllBytes(source))).append("\",\n");
        }
        json.append("  \"elapsed_seconds\": ").append((System.nanoTime() - started) / 1_000_000_000.0).append("\n");
        json.append("}\n");
        Files.writeString(output, json.toString());
        System.out.print(json);
        if (!uncovered.isEmpty()) {
            throw new IllegalStateException("Some targets had no allowed summand <=211");
        }
    }
}
