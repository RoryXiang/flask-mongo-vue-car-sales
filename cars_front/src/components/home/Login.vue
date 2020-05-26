<template>
  <div class="login_container">
      <login-top middTop="登录"></login-top>
      <div class="login_box">
        <el-form>
          <el-form-item>
            <el-input placeholder="帐号"></el-input>
          </el-form-item>
          <el-form-item>
            <el-input placeholder="密码"></el-input>
          </el-form-item>
          <el-form-item>
            <el-button type="primary">登录</el-button>
            <el-button type="info">重置</el-button>
          </el-form-item>
        </el-form>
      </div>
  </div>
</template>

<script>
import LoginTop from '@/components/commen/LoginTop.vue'
export default {
  components: {
      LoginTop
  }
}
</script>

<style>
.login_container {
    background-color: #2b4b6b;
    height: 100%;
    padding: 0px;
    border: 0px;
}
</style>