"""Model-coefficient helpers shared by main/supplementary plots."""

def get_task_slope_intercept(res, task, ref_task, measure):
    params = res.params
    inter_names = [
        f"{measure}:Task[T.{task}]",
        f"Task[T.{task}]:{measure}"
    ]

    beta0 = params.get("Intercept", 0.0)
    beta_measure = params.get(measure, 0.0)

    if task == ref_task:
        return beta0, beta_measure

    beta_task = params.get(f"Task[T.{task}]", 0.0)

    beta_int = 0.0
    for name in inter_names:
        if name in params.index:
            beta_int = params[name]
            break

    return beta0 + beta_task, beta_measure + beta_int
