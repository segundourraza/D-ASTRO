function [polyLow, polyUp] = computePolynomials(BC_range, gammalb, gammaub)
A = [ones([1,7]); BC_range;BC_range.^2; BC_range.^3; BC_range.^4; BC_range.^5; BC_range.^6; BC_range.^7]';
polyLow = flip((A\gammalb))';
polyUp = flip((A\gammaub))';