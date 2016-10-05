# Closest points in 3-D space
Course: Algorithms

Date: Fall 2015 (5th semester)

Find the pair of points that are closest to one another in 3-d space using a divide and conquer algorithm. The original assignment obfuscated this by describing it in terms of reactivity between different compounds, each made up of the elements N, C and O. The reactivity formula just happened to be the same as the formular for calculating the distance between two points in 3-d Euclidean space. Hence the class name: "BigBadaBoom"

Notes:
* We were required to submit a solution using a divide and conquer strategy.
* If I recall correctly, we needed to handle datasets of up to 10,000,000 points within 5 seconds.
* I didn't know it at the time, but this is one of the problems examined in Cormen's *Algorithms* book.
* It turns out I used the correct approach, but I needed to "merge" across one more dimension, so performance isn't quite as good as it could be. It's pretty close though--roughly O(n log(n))
* I would ordinarily refrain from writing 13 nested ternaries; in this case it was part of an inside joke with the professor teaching the course.

