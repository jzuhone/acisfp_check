from ..acisfp_check import model_path, ACISFPCheck
from acis_thermal_check.regression_testing import \
    RegressionTester, all_loads
import pytest


@pytest.fixture(autouse=True, scope='module')
def fp_rt(test_root):
    # ACIS state builder tests
    rt = RegressionTester(ACISFPCheck, model_path, "acisfp_test_spec.json",
                          test_root=test_root, sub_dir='acis')
    rt.run_models(state_builder='acis')
    return rt

# Prediction tests

@pytest.mark.parametrize('load', all_loads)
def test_prediction(fp_rt, answer_store, load):
    fp_rt.run_test("prediction", load, answer_store=answer_store)

# Validation tests

@pytest.mark.parametrize('load', all_loads)
def test_validation(fp_rt, answer_store, load):
    fp_rt.run_test("validation", load, answer_store=answer_store)
