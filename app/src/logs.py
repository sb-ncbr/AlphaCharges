class Logs:
    def __init__(self,
                 data_dir: str,
                 empirical_method: str):
        self.data_dir = data_dir
        self.step = 1
        self.max_steps = 6 if empirical_method == 'SQEqp' else 7

    def add_log(self,
                log: str):
        html_log = f'<p><strong> Step {self.step}/{self.max_steps}:</strong> {log}</p>\n'
        previous_logs = open(f'{self.data_dir}/page_log.txt').readlines()
        with open(f'{self.data_dir}/page_log.txt', 'w') as page_log_file:
            if len(previous_logs) and '...' in previous_logs[-1]:
                previous_logs = previous_logs[:-1]
                self.step += 1
            page_log_file.write(''.join(previous_logs) + html_log)