"""
User based Collaborative Filtering

dumb example

A Survey of Collaborative Filtering Techniques
Xiaoyuan Su and Taghi M. Khoshgoftaar
http://downloads.hindawi.com/archive/2009/421425.pdf

https://realpython.com/build-recommendation-engine-collaborative-filtering/
https://medium.com/@cfpinela/recommender-systems-user-based-and-item-based-collaborative-filtering-5d5f375a127f

Memory based? could use kmeans?
"""

using Distances;
using LinearAlgebra;

# number of items
N = 2

# number of users
U = 4

# utility matrix of users ratings of items

# so need a rating for user 1 item 2

user_item_matrix = [
        1 0 ; 
        2 4 ; 
        2.5 4 
        4.5 5 
    ]

# similar users: use either cosine or pearson corr
user_similarity = zeros(U,U)
for i in 1:U
    for j in 1:U
        user_similarity[i,j] = cosine_dist(user_item_matrix[i, :], user_item_matrix[j, :] )
    end
end 

# calculate rating user 1, item 2 
# using a wieghted average: sum(rating_item_i_by_user * similarity to target user) / sum(similarites of all users) 
rating = dot(user_similarity[:,1], user_item_matrix[ : , 2]) / sum(user_similarity[:,1])

@show rating
