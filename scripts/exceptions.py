"""自定义异常体系.

异常层次:
    SensoryGenomicsError
        ├── VcfError
        ├── VepError
        ├── ConfigError
        ├── ApiError
        └── AssessmentError
"""


class SensoryGenomicsError(Exception):
    """所有 Skill 相关异常的基类."""

    def __init__(self, message: str = "") -> None:
        super().__init__(message)
        self.message = message


class VcfError(SensoryGenomicsError):
    """VCF 文件解析或处理异常."""

    def __init__(self, message: str = "", record: str = "") -> None:
        super().__init__(message)
        self.record = record


class VepError(SensoryGenomicsError):
    """VEP 注释服务异常."""

    def __init__(self, message: str = "", status_code: int = 0) -> None:
        super().__init__(message)
        self.status_code = status_code


class ConfigError(SensoryGenomicsError):
    """配置加载或校验异常."""

    def __init__(self, message: str = "", field: str = "") -> None:
        super().__init__(message)
        self.field = field


class ApiError(SensoryGenomicsError):
    """外部 API 调用异常（UniProt/gnomAD/ClinVar/GTEx）."""

    def __init__(self, message: str = "", api_name: str = "", status_code: int = 0) -> None:
        super().__init__(message)
        self.api_name = api_name
        self.status_code = status_code


class AssessmentError(SensoryGenomicsError):
    """功能影响评估引擎异常."""

    def __init__(self, message: str = "", gene_symbol: str = "") -> None:
        super().__init__(message)
        self.gene_symbol = gene_symbol
