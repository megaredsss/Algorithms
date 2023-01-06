package com.mastermind.morgoon.AlgoSolver;

import java.util.ArrayList;
import java.util.Collections;


public class AlgoSolver {
    static class Complex {
        public final double re;
        public final double im;

        public Complex() {
            this(0, 0);
        }

        public Complex(double real, double imag) {
            re = real;
            im = imag;
        }

        public Complex plus(Complex b) {
            return new Complex(this.re + b.re, this.im + b.im);
        }

        public Complex minus(Complex b) {
            return new Complex(this.re - b.re, this.im - b.im);
        }

        public Complex multi(Complex b) {
            return new Complex(this.re * b.re - this.im * b.im,
                    this.re * b.im + this.im * b.re);
        }

        @Override
        public String toString() {
            return String.format("(%f,%f)", re, im);
        }
    public static int bitReverse(int n, int bits) {
        int reversedN = n;
        int count = bits - 1;

        n >>= 1;
        while (n > 0) {
            reversedN = (reversedN << 1) | (n & 1);
            count--;
            n >>= 1;
        }

        return ((reversedN << count) & ((1 << bits) - 1));
    }

    static void fft(Complex[] buffer) {

        int bits = (int) (Math.log(buffer.length) / Math.log(2));
        for (int j = 1; j < buffer.length / 2; j++) {

            int swapPos = bitReverse(j, bits);
            Complex temp = buffer[j];
            buffer[j] = buffer[swapPos];
            buffer[swapPos] = temp;
        }

        for (int N = 2; N <= buffer.length; N <<= 1) {
            for (int i = 0; i < buffer.length; i += N) {
                for (int j = 0; j < N / 2; j++) {

                    int evenIndex = i + j;
                    int oddIndex = i + j + (N / 2);
                    Complex EvenPart = buffer[evenIndex];
                    Complex OddPart = buffer[oddIndex];

                    double term = (-2 * Math.PI * j) / (double) N;
                    Complex exp = (new Complex(Math.cos(term), Math.sin(term)).multi(OddPart));

                    buffer[evenIndex] = EvenPart.plus(exp);
                    buffer[oddIndex] = EvenPart.minus(exp);
                }
            }
        }
    }
        public static int[][] FloydMarshall(int[][] graph) {
            final int INF = 999999;
            for (int k = 0; k < graph.length; k++) {
                for (int i = 0; i < graph.length; i++) {
                    for (int j = 0; j < graph.length; j++) {

                        if (graph[i][k] + graph[k][j] < graph[i][j])
                            graph[i][j] = graph[i][k] + graph[k][j];
                    }
                }
            }
            for (int i = 0; i < graph.length; ++i) {
                for (int j = 0; j < graph.length; ++j) {
                    if (graph[i][j] == INF)
                        System.out.print("INF ");
                    else
                        System.out.print(graph[i][j] + "  ");
                }
                System.out.println();
            }
            return graph;
        }

        public static ArrayList<Integer> increasingArrayList(ArrayList<Integer> numbers) {
            int max = 0;
            ArrayList<Integer> tmp = new ArrayList<>();
            for (int i = 0; i < numbers.size(); i++) {
                tmp.add(0);
            }

            ArrayList<Integer> result = new ArrayList<>();
            for (int i = 1; i < numbers.size(); i++) {
                if (numbers.get(i) > numbers.get(i - 1)) {
                    tmp.set(i, tmp.get(i - 1) + 1);
                } else tmp.set(i, 0);
                {
                    max = Math.max(tmp.get(i), max);
                }
            }
            for (int i = 0; i < numbers.size(); i++) {
                if (tmp.get(i) == max) {
                    for (int j = i - tmp.get(i); j <= i; j++) {
                        result.add(numbers.get(j));
                    }
                }
            }
            return result;
        }

        public static ArrayList<Integer> decreasingArrayList(ArrayList<Integer> numbers) {
            int max = 0;
            ArrayList<Integer> tmp = new ArrayList<>();
            for (int i = 0; i < numbers.size(); i++) {
                tmp.add(0);
            }
            ArrayList<Integer> result = new ArrayList<>();
            for (int i = 1; i < numbers.size(); i++) {
                if (numbers.get(i) < numbers.get(i - 1)) {
                    tmp.set(i, tmp.get(i - 1) + 1);
                } else tmp.set(i, 0);
                {
                    max = Math.max(tmp.get(i), max);
                }
            }

            if (max == 0) {
                result = null;
            } else {
                for (int i = 0; i < numbers.size(); i++) {
                    if (tmp.get(i) == max) {
                        for (int j = i - tmp.get(i); j <= i; j++) {
                            result.add(numbers.get(j));
                        }
                    }
                }
            }
            return result;
        }

        public static float median(ArrayList<Integer> numbers) {
            float median = 0;
            Collections.sort(numbers);
            if (numbers.size() % 2 == 0) {
                median = ((float) numbers.get(Math.round(((float) numbers.size()) / 2)) + ((float) numbers.get(Math.round(((float) numbers.size()) / 2) - 1))) / 2;
            } else {
                median = numbers.get(numbers.size() / 2);


            }
            return median;
        }

        public static void main(String[] args) {
            ArrayList<Integer> testArrayList = new ArrayList<>();
            for (int i = 0; i < 10; i++) {
                testArrayList.add(i);
            }

            ArrayList<Integer> result = increasingArrayList(testArrayList);
            System.out.println(result);
            result = decreasingArrayList(testArrayList);
            System.out.println(result);
            float median = median(testArrayList);
            System.out.println(median);
            int[][] graph = {{0, 5, 999999, 10}, {999999, 0, 3, 999999}, {999999, 999999, 0, 1}, {999999, 999999, 999999, 0}};
            graph = FloydMarshall(graph);
            double[] signal = {9.0, 1.0, 4.0, 2.6, 81.6, 13.9, 25.0, 0.0};

            Complex[] testsignal = new Complex[signal.length];
            for (int i = 0; i < signal.length; i++)
                testsignal[i] = new Complex(signal[i], 0.0);

            fft(testsignal);
            for (Complex resultFFt : testsignal) {
                System.out.println(resultFFt);
            }

        }
    }
}