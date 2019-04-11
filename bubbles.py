import random

goods = 20

# number of goods to select in the simulation?
number_of_selections = 10


# TODO how to get actual ex_ante probability
current_beliefs = [random.uniform(0,1) for i in range(goods)]

# To keep track of what goods we chosen so far
choices = [i+1 for i in range(goods)]


def make_consumer_choice(current_beliefs):
    """ 
    """
    # choose a good
    # TODO how to choose
    index = random.randint(0,len(choices)) - 1
    # pop chosen good dont want to choose this good again
    choice = choices.pop(index)
    print("chose good #: {}".format(choice))
    #print(choices)
    
    #TODO an actual update current_beliefs func
    map(lambda x : x-.00001, current_beliefs) 

# go through select unique items, update beliefs 
for i in range(number_of_selections):
    print("iteration: {}".format(i))
    print(current_beliefs)
    make_consumer_choice(current_beliefs)

print("final beliefs: {}".format(current_beliefs))
