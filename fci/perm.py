#!/usr/bin/python
'''From here: http://stackoverflow.com/questions/4250125/generate-permutations-of-list-with-repeated-elements
and here: http://blog.bjrn.se/2008/04/lexicographic-permutations-using.html'''
 
def lexicographic_permutations(seq):
    '''Algorithm L for The Art of Computer Programming, Volume 4, Fascicle 2: Generating All Tuples and Permutations.'''
 
    def reverse(seq, start, end):
        '''In-place reversion. Start and end denote the slice ends which forms the chunk to reverse'''
        # seq = seq[:start] + reversed(seq[start:end]) + \
        #       seq[end:]
        
        end -= 1
        while start < end:
            seq[start], seq[end] = seq[end], seq[start]
            start += 1
            end -= 1
 
 
    #Some checks
    if not seq:
        raise StopIteration
 
    try:
        seq[0]
    except TypeError:
        raise TypeError("seq must allow random access.")
 
 
    end = len(seq)
    seq = seq[:] #copy the input sequence
    
    yield seq #return the seq itself as a first value
 
   
    while True:
        j = end - 1
        while True:
            # Step 1.
            j -= 1
            if j == -1:
                raise StopIteration 
            if seq[j] < seq[j+1]:
                # Step 2.
                l = end - 1
                while not seq[j] < seq[l]:
                    l -= 1
                seq[j], seq[l] = seq[l], seq[j]
                # Step 3.
                reverse(seq, j+1, end)
                # Change to yield references to get rid of
                # (at worst) |seq|! copy operations.
                yield seq[:]
                break
 
 
if __name__ == '__main__':
    import timeit
    #use 'pass' instead of 'print(p)' to get rid of hefty print function time complexity charges.
    t = timeit.Timer(stmt='for p in lexicographic_permutations(range(9)): pass',setup='from __main__ import lexicographic_permutations')
    print(t.timeit(number=1))

    for p in lexicographic_permutations([1, -2, -2]):
        print p