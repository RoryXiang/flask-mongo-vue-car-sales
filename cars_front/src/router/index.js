import Vue from 'vue'
import Router from 'vue-router'
import Login from '@/components/home/Login.vue'
import Register from '@/components/admin/Register.vue'

Vue.use(Router)

export default new Router({
  routes: [
    {
      path: '/login',
      component: Login 
    },
    {
      path: '/register',
      component: Register
    }
  ]
})
