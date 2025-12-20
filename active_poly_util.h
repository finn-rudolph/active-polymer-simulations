inline double correct_coord_diff(double diff, double box_len) {
    if (diff > box_len / 2)
        return diff - box_len;
    else if (diff < -box_len / 2)
        return diff + box_len;
    return diff;
}

inline double correct_coord(double coord, double reference, double box_len) {
    double difference = coord - reference;
    if (difference > box_len / 2)
        return coord - box_len;
    else if (difference < -box_len / 2)
        return coord + box_len;
    return coord;
}