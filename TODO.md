## TODO

A few things need to happen for partitioning euclidean graphs. Not all of these
may be needed, but these are the places I should look into.

- [ ] Edge matching
    - matching phase of coarsening is done using edge weights
    - I'm not sure how to best include location in 3D space
    - Probably the lowest priority

- [ ] Initial partitioning
    - Probably not much needed for this step, since it just invokes recursive
      bisectioning

- [ ] Uncoarsening
    - Probably the place to look
    - Works off of two priority queues and competing gain + locking/unlocking of
      vertices
    - Still going to be challenging to work Euclidean metrics into priority
      queue
    - Maybe there's a paper on that
