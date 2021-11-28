#pragma once

inline static double weight(double distance, double re){return (distance >= re)?(0):((re / distance) - 1.0);}
