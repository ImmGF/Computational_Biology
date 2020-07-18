Implement the simplified consensus method for motif finding (based on slide 20) and use it on the E.coli  promoter sequences 
for motifs of length 3,4,5,6. The algorithm works as follows:

- based on any position in the first promoter sequence initialize a PWM
- for every following step:

      a) choose the position  in one of the remaining sequences that is best matching the current PWM

      b) modify your PWM, by adding the word at the chosen position to your PWM and remove this sequence from further considerations

- Perform this procedure for all possible starting points in the first sequence and return a few  (chosen by the user, by default 5) resulting PWMs that have the highest information content
    
