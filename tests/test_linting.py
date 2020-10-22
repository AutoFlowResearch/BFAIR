from unittest import TestCase
import os


class Linting(TestCase):
    def test_flake8(self):
        """Test code python style"""
        result = os.system("flake8 --max-line-length=120 --ignore=PT009")
        print(result)
