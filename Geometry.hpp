class EucludianSpace
{
    /*
    * The CartesianPoint structure is only intended to be
    * used when there's a need for exact cartesian coordinates
    * of a certain object (a node, a drone or anything else)
    * 
    * It is implied that most of the computations considering
    * topology and euclidian distance are implemented internally
    * within the classes that provide particular mechanisms for
    * their Euclidian space and DO NOT REQUIRE their internal
    * coordinates (that, by the way, can be integers) to be converted
    * to cartesian.
    */
    template<unsigned int Dimensions>
    struct CartesianPoint
    {
        double Coordinates[Dimensions];
    };

    class IGrid
    {

    };
};