def distance_func(x, y):
    return x / y


class FoldChangeCalculator:
    def __init__(self, leftCondition, rightCondition):
        self.leftCondition = leftCondition
        self.rightCondition = rightCondition

    def calculate_fold_change(self):
        leftAverage = self.leftCondition.get_averages_over_repetitions()
        rightAverage = self.rightCondition.get_averages_over_repetitions()
        return distance_func(leftAverage, rightAverage)

